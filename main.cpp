#include <chrono>
#include <fstream>

#include <spdlog/spdlog.h>
#include <htslib/vcf.h>
#include <htslib/hts.h>

#include "Params.hpp"
#include "GlobalStore.hpp"
#include "processLine.hpp"
#include "paramsControl.hpp"
#include "threadPool.hpp"


std::queue<bcf1_t*> writeQueue;
std::mutex queueMutex;
std::condition_variable cv;
bool done = false;
// 写入线程
void writerThread(htsFile* out_file, bcf_hdr_t* header) {
    while (true) {
        std::queue<bcf1_t *>::value_type record;
        {
            std::unique_lock lock(queueMutex);
            cv.wait(lock, [] { return !writeQueue.empty() || done; });
            if (done && writeQueue.empty()) break;
            record = writeQueue.front();
            writeQueue.pop();
        }
        if (record) {
            if (bcf_write(out_file, header, record) != 0) {
                spdlog::error("Failed to write record");
            }
            bcf_destroy(record);
        }
    }
}

int main(const int argc, const char *argv[]) {
    // 记录程序开始时间
    const auto start = std::chrono::high_resolution_clock::now();
    spdlog::set_level(spdlog::level::info);
    spdlog::info("Start!");
    Params params;
    Controller controller{};
    checkParameters(argc, argv, params);
    // 打开 VCF、BCF 文件
    htsFile *input_file = hts_open(params.inputFilename.c_str(), "r");
    if (!input_file) {
        spdlog::error("Failed to open input file:{}", params.inputFilename);
        return 1;
    }

    // 打开输出文件
    htsFile *output_file = hts_open(params.outputFilename.c_str(), "wz");
    if (!output_file) {
        spdlog::error("Unable to create file: {}, please check permissions", params.outputFilename);
        return 1;
    }

    // 读取 VCF 文件头部
    bcf_hdr_t *header = bcf_hdr_read(input_file);

    // 启动写入线程
    std::thread writer(writerThread, output_file, header);

    // 获取样本数量
    controller.n_samples = bcf_hdr_nsamples(header);
    spdlog::info("N sample: {}", controller.n_samples);

    bcf_hdr_remove(header, BCF_HL_FMT, "GQ");
    bcf_hdr_remove(header, BCF_HL_FMT, "PGT");
    bcf_hdr_remove(header, BCF_HL_FMT, "PID");
    bcf_hdr_remove(header, BCF_HL_FMT, "PL");
    bcf_hdr_remove(header, BCF_HL_FMT, "PS");
    bcf_hdr_remove(header, BCF_HL_FMT, "MIN_DP");
    bcf_hdr_remove(header, BCF_HL_FMT, "RGQ");
    bcf_hdr_remove(header, BCF_HL_FMT, "SB");
    bcf_hdr_remove(header, BCF_HL_FLT, nullptr);

    for (int i = 0; i < header->n[BCF_DT_ID]; ++i) {
        const char* key = bcf_hdr_int2id(header, BCF_DT_ID, i);
        // 只保留指定的 INFO 字段
        if (bcf_hdr_idinfo_exists(header, BCF_HL_INFO, i)) {
            if (!params.keep_fields.contains(key)) {
                bcf_hdr_remove(header, BCF_HL_INFO, key);
            }
        }
    }

    // 写入头部到输出文件
    if (bcf_hdr_write(output_file, header) < 0) {
        bcf_hdr_destroy(header);
        hts_close(input_file);
        hts_close(output_file);
        spdlog::error("Can not write header to file {}!", params.outputFilename);
        return 1;
    }

    // 创建记录对象
    bcf1_t *record = bcf_init();
    if (!record) {
        spdlog::error(" Could not initialize record.");
        bcf_hdr_destroy(header);
        bcf_close(input_file);
        return 1;
    }

    // 初始化计数器
    int record_counter = 0;

    // 逐行读取 VCF 文件
    while (bcf_read(input_file, header, record) == 0) {
        // 解析等位基因数
        controller.n_alleles = record->n_allele;
        bcf1_t *local_record = bcf_dup(record);
        handleRecord(header, local_record, params, controller, std::ref(writeQueue), std::ref(queueMutex), std::ref(cv));
        // 每处理 20000 行，输出当前染色体编号和位置
        if (++record_counter % 20000 == 0) {
            spdlog::info("Current location: {}:{}", bcf_hdr_id2name(header, record->rid), record->pos + 1);
        }
    }

    // 结束写入线程
    {
        std::lock_guard lock(queueMutex);
        done = true;
    }
    cv.notify_one();
    writer.join();

    // 清理资源
    bcf_hdr_destroy(header);
    hts_close(input_file);
    hts_close(output_file);

    // **创建 VCF.gz 索引**
    if (bcf_index_build(params.outputFilename.c_str(), 0) < 0) {
        spdlog::error("Failed to create index!");
        return 1;
    }

    // 记录程序结束时间
    const auto end = std::chrono::high_resolution_clock::now();
    // 计算时间差（以毫秒为单位）
    const auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    spdlog::info("Total running time: {:.3f}s", duration.count());
    return 0;
}