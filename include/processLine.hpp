//
// Created by 王亚洲 on 25-1-2.
//

#ifndef PROCESSLINE_HPP
#define PROCESSLINE_HPP

#include "Params.hpp"
#include "GlobalStore.hpp"
#include <htslib/vcf.h>
#include <threadPool.hpp>

void handleRecord(const bcf_hdr_t *header, bcf1_t *record, const Params &params, const Controller &controller, std::queue<bcf1_t *> &writeQueue, std
                  ::mutex &queueMutex, std::condition_variable &cv);

#endif //PROCESSLINE_HPP
