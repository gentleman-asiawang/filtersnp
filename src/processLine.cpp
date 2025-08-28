//
// Created by 王亚洲 on 25-1-2.
//
#include "processLine.hpp"
#include "RewriteRecord.hpp"
#include "spdlog/spdlog.h"

#include <algorithm>

void handleRecord(const bcf_hdr_t *header, bcf1_t *record, const Params &params, const Controller &controller, std::queue<bcf1_t *> &writeQueue, std
                  ::mutex &queueMutex, std::condition_variable &cv) {
    spdlog::debug("Record: {}:{}, n_alleles:{}", bcf_hdr_id2name(header, record->rid), record->pos + 1, controller.n_alleles);
    // 根据QUAL过滤
    if (record->qual < 30.0) return;
    bcf_unpack(record, BCF_UN_INFO);
    //  根据QD字段过滤
    float* qd_value = nullptr;
    int nqd = 0;
    bcf_get_info_float(header, record, "QD", &qd_value, &nqd);

    if (nqd >0 && qd_value[0] < 2.0) return;
    //  根据FS字段过滤
    float* fs_value = nullptr;
    int nfs = 0;
    bcf_get_info_float(header, record, "FS", &fs_value, &nfs);
    float* rprs_value = nullptr;
    int nrprs = 0;
    // 根据ReadPosRankSum字段过滤
    bcf_get_info_float(header, record, "ReadPosRankSum", &rprs_value, &nrprs);
    if (params.filetype == "SNP") { // SNP
        if (nfs > 0 && fs_value[0] > 60.0) return;
        if (nrprs > 0 && rprs_value[0] < -8.0) return;
        // 根据MQ字段过滤
        float* mq_value = nullptr;
        int nmq = 0;
        bcf_get_info_float(header, record, "MQ", &mq_value, &nmq);
        if (nmq > 0 && mq_value[0] < 40.0) return;
        // 根据MQRankSum字段过滤
        float* MQRankSum_value = nullptr;
        int nMQRankSum = 0;
        bcf_get_info_float(header, record, "MQRankSum", &MQRankSum_value, &nMQRankSum);
        if (nMQRankSum > 0 && MQRankSum_value[0] < -12.5) return;
    } else {
        if (nfs > 0 && fs_value[0] > 200.0) return;
        if (nrprs>0 && rprs_value[0] < -20.0) return;
    }

    // 提取 AF (Allele Frequency) 字段
    float* af_values = nullptr;
    int naf = 0; // 初始化数组大小
    bcf_get_info_float(header, record, "AF", &af_values, &naf);
    if (naf > 0) {
        // 使用 std::all_of 检查是否所有 af_values 的值都小于 0.01
        if (std::all_of(af_values, af_values + naf, [](const double af) { return af <= 0.01; })) {
            return;  // 如果所有值都小于等于 0.01，不输出这行
        }
        // 着一行不因AF删除，则先计算缺失率，根据缺失率判断
        // 获取AN字段
        int *an_value = nullptr;
        int nan = 0;
        bcf_get_info_int32(header, record, "AN", &an_value, &nan);
        if (const float missing_rate = 1 - static_cast<float>(an_value[0]) / static_cast<float>(controller.n_samples) * 2; missing_rate > 0.01) return;
    } else {
        spdlog::warn("AF field not found or empty.");
        return;
    }

    // record->n_info = 0;    // 更新记录中 INFO 字段为NULL

    // 遍历每条记录
    // 遍历所有 INFO 字段
    for (int i = 0; i < header->n[BCF_DT_ID]; ++i) {
        const char* key = bcf_hdr_int2id(header, BCF_DT_ID, i);

        // 检查字段是否为 INFO 类型
        if (bcf_hdr_idinfo_exists(header, BCF_HL_INFO, i)) {
            // 如果不是需要保留的字段，删除它
            if (!params.keep_fields.contains(key)) {
                bcf_update_info(header, record, key, nullptr, 0, 0);
            }
        }
    }

    // 删除一些用不到的信息
    if (bcf_update_format_float(header, record, "GQ", nullptr, 0) != 0) {
        spdlog::error("Failed to clear format field");
    }
    if (bcf_update_format_float(header, record, "PGT", nullptr, 0) != 0) {
        spdlog::error("Failed to clear format field");
    }
    if (bcf_update_format_float(header, record, "PID", nullptr, 0) != 0) {
        spdlog::error("Failed to clear format field");
    }
    if (bcf_update_format_float(header, record, "PL", nullptr, 0) != 0) {
        spdlog::error("Failed to clear format field");
    }
    if (bcf_update_format_float(header, record, "MIN_DP", nullptr, 0) != 0) {
        spdlog::error("Failed to clear format field");
    }
    if (bcf_update_format_float(header, record, "RGQ", nullptr, 0) != 0) {
        spdlog::error("Failed to clear format field");
    }
    if (bcf_update_format_float(header, record, "SB", nullptr, 0) != 0) {
        spdlog::error("Failed to clear format field");
    }
    if (bcf_update_format_float(header, record, "PS", nullptr, 0) != 0) {
        spdlog::error("Failed to clear format field");
    }
    {
        // 加锁并加入队列
        std::lock_guard<std::mutex> lock(queueMutex);
        writeQueue.push(record);
    }
    cv.notify_one();  // 通知写入线程
}