# rm Parallel_CaLiG

icpx parallel_calig.cpp -o Parallel_CaLiG -O3 -qopenmp -std=c++11

# ./Parallel_CaLiG -d /home/cc/haibin2/CSM-Benchmark/ContinuousSubgraphMatching/matching/CaLiG/CaLiG/github10/initial \
# -s /home/cc/haibin2/CSM-Benchmark/ContinuousSubgraphMatching/matching/CaLiG/CaLiG/github10/s \
# -q CSM-Benchmark/ContinuousSubgraphMatching/matching/CaLiG/CaLiG/github10/Q/6/q2 -t 1

 ./Parallel_CaLiG -d /home/cc/haibin2/CSM-Benchmark/ContinuousSubgraphMatching/matching/CaLiG/CaLiG/github10/initial \
 -s /home/cc/haibin2/CSM-Benchmark/ContinuousSubgraphMatching/matching/CaLiG/CaLiG/github10/s \
  -q /home/cc/haibin2/CSM-Benchmark/ContinuousSubgraphMatching/matching/CaLiG/CaLiG/github10/Q/6/q1 -t 4


./Parallel_CaLiG -d /home/cc/haibin2/livejournal/data_graph/data_Unlabel.graph \
-s /home/cc/haibin2/livejournal/data_graph/insertion_Unlabel.graph \
-q /home/cc/haibin2/livejournal/unlabel_query/8_self/dense/Q_11