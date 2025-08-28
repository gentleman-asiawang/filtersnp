// Created by 王亚洲 on 25-1-2.
//

#ifndef REWRITERECORD_HPP
#define REWRITERECORD_HPP

#include <htslib/vcf.h>
#include "GlobalStore.hpp"

void RewriteRecord(const bcf_hdr_t *header, bcf1_t *record, int32_t *gt_arr, int ngt_arr, const Controller& controller, int drop_idx);

#endif //REWRITERECORD_HPP
