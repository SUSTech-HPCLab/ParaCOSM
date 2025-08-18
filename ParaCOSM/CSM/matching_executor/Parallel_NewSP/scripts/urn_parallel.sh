#!/bin/bash
# 批量测试从 Q_1 到 Q_100 的所有图

export LD_LIBRARY_PATH=/home/cc/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

# 数据目录
DIR=/home/cc/haibin2/
# 数据集
DATA_SET=amazon/6
# 时间限制
TIME_LIMIT=3600

# 创建结果目录
RESULT_DIR="results_$(date +%Y%m%d_%H%M%S)"
mkdir -p $RESULT_DIR

# 数据图和插入图路径
DATA_GRAPH="/home/cc/haibin2/amazon/6/data_graph/data.graph"
INSERT_GRAPH="/home/cc/haibin2/amazon/6/data_graph/insertion.graph"
QUERY_BASE_DIR="/home/cc/haibin2/amazon/8_self/sparse"

echo "开始批量测试 Q_1 到 Q_100..."
echo "结果将保存在 $RESULT_DIR 目录下"

# 遍历从 Q_1 到 Q_100
for i in $(seq 1 100); do
    QUERY_FILE="${QUERY_BASE_DIR}/Q_${i}"
    
    # 检查文件是否存在
    if [ -f "$QUERY_FILE" ]; then
        echo "测试 Q_${i}..."
        
        # 运行命令并将输出保存到文件
        OUTPUT_FILE="${RESULT_DIR}/result_Q_${i}.txt"
        sudo ./build/csm -d $DATA_GRAPH -u $INSERT_GRAPH \
             -q $QUERY_FILE --time-limit ${TIME_LIMIT} --report-initial off > $OUTPUT_FILE 2>&1
        
        # 提取关键结果信息
        echo "Q_${i} 测试完成。结果保存在 $OUTPUT_FILE"
    else
        echo "警告: $QUERY_FILE 不存在，跳过。"
    fi
done

echo "所有测试完成！"
echo "结果位于: $RESULT_DIR"