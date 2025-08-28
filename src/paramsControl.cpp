//
// Created by 王亚洲 on 24-1-11.
//

#include "paramsControl.hpp"
#include <iostream>
#include <cstdlib>
#include <spdlog/spdlog.h>

static void printUsage(const std::string &programName) {
    std::cout << "Usage: " << programName << " -i example.vcf.gz -o example.vcf.gz" << std::endl;
    std::cout << "This program is used to filter vcf files" << std::endl; // 描述用途
    std::cout << "Options:" << std::endl;
    std::cout << "  -i                  <filepath>      Required: Input file (support [vcf, bcf, vcf.gz, bcf.gz], Automatically determine based on file suffix)." << std::endl;
    std::cout << "  -o                  <filepath>      Required: Output file (support [vcf, bcf, vcf.gz, bcf.gz], Automatically determine based on file suffix)." << std::endl;
    std::cout << "  -type               <SNP;INDEL>     Required: One of {SNP, INDEL}." << std::endl;
    std::cout << "  --version                           Show this version message" << std::endl;
    std::cout << "  -h                                  Show this help message" << std::endl;
}

void printVersion() {
    std::cout << "Filtersnp\nVersion: V-1.0.0\tLast update:2024.1.13" << std::endl;
}

void checkParameters(const int& argc, const char *argv[], Params& params) {
    // 检查是否提供了任何参数
    if (argc == 1) { // 只有程序名称，没有其他参数
        spdlog::error("No parameters provided. Use -h for help.");
        printUsage(argv[0]); // 打印帮助信息
        std::exit(1);        // 退出程序
    }

    // 解析参数，支持任意顺序
    for (int i = 1; i < argc; ++i) {
        if (std::string arg = argv[i]; arg == "-i") {
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                params.inputFilename = argv[++i]; // 获取输入文件名
            } else {
                spdlog::error("No input filename specified after -i");
                std::exit(1);
            }
        } else if (arg == "-o") {
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                params.outputFilename = argv[++i]; // 获取输出文件名
            } else {
                spdlog::error("No output filename specified after -o");
                std::exit(1);
            }
        }else if (arg == "-type") {
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                if (std::string(argv[i+1]) != "SNP" && std::string(argv[i+1]) != "INDEL") {
                    spdlog::error("Please input SNP or INDEL");
                    std::exit(1);
                }
                params.filetype = argv[++i];
            } else {
                spdlog::error("Nothing specified after -type");
                std::exit(1);
            }
        } else if (arg == "-h") {
            printUsage(argv[0]); // 打印使用说明
            std::exit(1);
        } else if (arg == "--version") {
            printVersion();
            std::exit(1);
        }
    }

    // 检查必须的参数是否已提供
    if (params.inputFilename.empty()) {
        spdlog::error("Missing required parameter: -i (input filename)");
        printUsage(argv[0]);
        std::exit(1);
    }
    if (params.outputFilename.empty()) {
        spdlog::error("Missing required parameter: -o (output filename)");
        printUsage(argv[0]);
        std::exit(1);
    }
    if (params.filetype.empty()) {
        spdlog::error("Missing required parameter: -type (SNP、INDEL)");
        printUsage(argv[0]);
        std::exit(1);
    }
}
