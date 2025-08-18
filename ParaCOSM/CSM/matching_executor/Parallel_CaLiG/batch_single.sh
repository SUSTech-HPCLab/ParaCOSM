#!/bin/bash
# filepath: /home/cc/haibin2/test_all_8_self_queries.sh

# 数据图和插入图路径
DATA_GRAPH="/home/cc/haibin2/livejournal/data_graph/data_Unlabel.graph"
INSERTION_GRAPH="/home/cc/haibin2/livejournal/data_graph/insertion_Unlabel.graph"


# export OMP_PROC_BIND=close
# export OMP_PLACES=threads
# export OMP_WAIT_POLICY=ACTIVE

# 查询目录的基础路径
QUERY_BASE_DIR="/home/cc/haibin2/livejournal/unlabel_query/8_self"

# 确保Parallel_CaLiG已经编译
cd /home/cc/haibin2/CSM-Benchmark/ContinuousSubgraphMatching/matching/CaLiG/
# icpx parallel_calig.cpp -o Parallel_CaLiG -O2 -qopenmp -std=c++11

# 创建结果目录
RESULTS_DIR="/home/cc/haibin2/8_self_single_results"
mkdir -p "$RESULTS_DIR"

# 记录开始时间
echo "开始测试时间: $(date)" > "$RESULTS_DIR/summary.txt"

# 遍历8_self下所有子目录
for SUB_DIR in "$QUERY_BASE_DIR"/*; do
    if [ -d "$SUB_DIR" ]; then
        SUB_DIR_NAME=$(basename "$SUB_DIR")
        echo "处理目录: $SUB_DIR_NAME"
        mkdir -p "$RESULTS_DIR/$SUB_DIR_NAME"
        
        # 遍历每个子目录中的所有Q文件
        for QUERY_FILE in "$SUB_DIR"/Q_*; do
            if [ -f "$QUERY_FILE" ]; then
                QUERY_NAME=$(basename "$QUERY_FILE")
                echo "测试查询: $QUERY_NAME"
                
                # 运行Parallel_CaLiG并记录输出
                echo "运行: ./CaliG_ori -d $DATA_GRAPH -s $INSERTION_GRAPH -q $QUERY_FILE"
                
                # 记录开始时间
                START_TIME=$(date +%s.%N)
                
                # 执行命令并捕获输出
                timeout 3600 ./CaliG_ori -d "$DATA_GRAPH" -s "$INSERTION_GRAPH" -q "$QUERY_FILE" > "$RESULTS_DIR/$SUB_DIR_NAME/$QUERY_NAME.log" 2>&1
                
                # 记录结束时间和耗时
                END_TIME=$(date +%s.%N)
                ELAPSED_TIME=$(echo "$END_TIME - $START_TIME" | bc)
                
                # echo "完成: $QUERY_NAME, 耗时: $ELAPSED_TIME 秒"
                echo "$SUB_DIR_NAME/$QUERY_NAME: $ELAPSED_TIME 秒" >> "$RESULTS_DIR/summary.txt"
            fi
        done
    fi
done

echo "所有测试完成: $(date)" >> "$RESULTS_DIR/summary.txt"
echo "结果保存在: $RESULTS_DIR"