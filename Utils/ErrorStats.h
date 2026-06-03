#pragma once

struct ErrorStats
{
    float maxFlat   = 0;
    float maxNagata = 0;
    double avgFlat   = 0;   
    double avgNagata = 0;
    long long samples = 0;
    float ratio = 0;        
};