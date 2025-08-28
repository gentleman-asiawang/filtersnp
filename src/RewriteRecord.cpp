//
// Created by 王亚洲 on 25-1-2.
//
#include "RewriteRecord.hpp"
#include <spdlog/spdlog.h>


void RewriteRecord(const bcf_hdr_t *header, bcf1_t *record, int32_t *gt_arr, int ngt_arr, const Controller& controller, const int drop_idx) {
    // spdlog::debug("RewriteRecord: {}:{}, n_alleles:{}, drop_idx:{}", bcf_hdr_id2name(header, record->rid), record->pos + 1, controller.n_alleles, drop_idx);
    // spdlog::debug("drop_alleles:{}", record->d.allele[drop_idx]);
    // spdlog::info("Record: {}:{}, n_alleles:{}, drop_idx:{}, drop_alleles:{}", bcf_hdr_id2name(header, record->rid), record->pos + 1, controller.n_alleles, drop_idx, record->d.allele[drop_idx]);
    // 去掉指定索引的alt
    if (drop_idx < 0 || drop_idx >= controller.n_alleles) {
        spdlog::error("dropIdx {:d} out of range", drop_idx);
        std::exit(EXIT_FAILURE); // 非零值表示程序异常退出
    }
    // 将后续的元素向前移动
    for (int i = drop_idx; i < record->n_allele - 1; i++) {
        record->d.allele[i] = record->d.allele[i + 1];
    }

    // 释放最后一个 allele
    record->d.allele[record->n_allele - 1] = nullptr;

    // 减少 n_allele 的数量
    record->n_allele--;

    // 获取GT
    bcf_get_format_int32(header, record, "GT", &gt_arr, &ngt_arr);
    for(int i=0; i < controller.n_samples; i++) {
        if (bcf_gt_allele(gt_arr[i * 2]) == drop_idx || bcf_gt_allele(gt_arr[i * 2 + 1]) == drop_idx) {
            // spdlog::debug("Delete idx: {}", drop_idx);
            // spdlog::debug("the:{}, gt_arr: {:#08b}/{:#08b}, GT: {}/{}",i , gt_arr[i * 2], gt_arr[i * 2 +1],  bcf_gt_allele(gt_arr[i * 2]), bcf_gt_allele(gt_arr[i * 2 + 1]));
            gt_arr[i * 2] = 0;
            gt_arr[i * 2 + 1] = 0;
        } else {
            if (bcf_gt_allele(gt_arr[i * 2]) > drop_idx) {
                gt_arr[i * 2] = bcf_gt_unphased(bcf_gt_allele(gt_arr[i * 2]) - 1);
            }
            if (bcf_gt_allele(gt_arr[i * 2 + 1]) > drop_idx) {
                gt_arr[i * 2 + 1] = bcf_gt_unphased(bcf_gt_allele(gt_arr[i * 2 + 1]) - 1);
            }
        }
    }


    // 更新AD
    // int *ad_arr = nullptr, ad_len = 0;
    // bcf_get_format_int32(header, record, "AD", &ad_arr, &ad_len);
    // spdlog::debug("ad_len: {:d}", ad_len);

    // // 创建新的 AD 数组
    // std::vector new_ad(n_samples * 2, 0);
    // if (!dropDP) {
    //     for(int i=0; i<n_samples; i++) {
    //         int dp = 0;
    //         for(int j=0; j<n_alleles; j++) {
    //             spdlog::debug("Nallele: {}, Nad: {}, Ngl: {}, sample: {}, ad{}: {}", n_alleles, ad_len, gl_len , i, j, ad_arr[i*n_alleles+j]);
    //             if (j != keep_idx1) {
    //                 if (ad_arr[i*n_alleles+j] > 0) dp += ad_arr[i*n_alleles+j];
    //             } else {
    //                 new_ad[i*2] = ad_arr[i*n_alleles+j];
    //             }
    //         }
    //         new_ad[i*2+1] = dp;
    //     }
    // } else {
    //     for(int i=0; i<n_samples; i++) {
    //         for(int j=0; j<n_alleles; j++) {
    //             if (j == keep_idx1) {
    //                 new_ad[i*2] = ad_arr[i*n_alleles+j];
    //             } else if (j == keep_idx2) {
    //                 new_ad[i*2+1] = ad_arr[i*n_alleles+j];
    //             }
    //         }
    //     }
    // }
    // // 更新记录中的 AD 数据
    // if (bcf_update_format_int32(header, record, "AD", new_ad.data(), new_ad.size()) != 0) {
    //     spdlog::error("Failed to update AD field");
    // }
    // free(ad_arr);
}