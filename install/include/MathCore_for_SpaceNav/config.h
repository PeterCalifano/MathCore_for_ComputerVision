/**
 * @file config.h.in
 * @author PeterC (petercalifano.gs@gmail.com)
 * @brief Configuration header file for library src, with variables set by CMake
 * @version 0.1
 * @date 2025-01-04
 */
#pragma once
#include <cstdio>
#include <string>

// Library version
#define PROJECT_VERSION_MAJOR 0
#define PROJECT_VERSION_MINOR 1
#define PROJECT_VERSION_PATCH 0
#define PROJECT_VERSION "0.1.0"
#define FULL_VERSION "0.1.0+4fe6965"


void PrintVersion()
{   
    printf("Application Version: %s\n", PROJECT_VERSION);
    printf("Full Version: %s\n", FULL_VERSION);
    printf("Version Details:\n");
    printf("Major: %d\n", PROJECT_VERSION_MAJOR);
    printf("Minor: %d\n", PROJECT_VERSION_MINOR);
    printf("Patch: %d\n", PROJECT_VERSION_PATCH);
    printf("Build Date: %s\n", __DATE__);
    printf("Build Time: %s\n", __TIME__);
}

std::string GetVersionString()
{
    return std::string(PROJECT_VERSION);
}

/* #undef VARIABLE_DEFINED_BY_CMAKE */

