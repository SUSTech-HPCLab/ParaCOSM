#include<iostream>
#include "core/StreamExecutionEnvironment.h"
#include "core/StreamInfo.h"
#include "core/DataGraph.h"
#include "core/BacktrackAlgorithm.h"
#include "core/DataStream.h"
#include "core/UpdateGraphSource.h"
#include "core/GraphSource.h"
#include "core/Algorithm.h"

// #include "ParaCSM.h"

int main(int argc, char** argv)
{

    StreamExecutionEnvironment env = StreamExecutionEnvironment.getExecutionEnvironment();

    StreamInfo streamInfo = env.loadStreamInfo("streamInfo.json");

    DataGraph graph = env.addDataGraphSource(streamInfo.getGraphSource());

    DataStream result = env.addUpdateGraphSource(streamInfo.getUpdateGraphSource());

    BacktrackAlgorithm algo = env.addBacktrackAlgorithm(streamInfo.getAlgorithm());

    auto result = env.execute(streamInfo.name());

    std::cout << "Result: " << result << std::endl;

    return 0;
}