#!/bin/bash
echo "----start shell----"
echo "开始编译generate_data.c文件"
gcc generate_data.c -o generate_data -g -Wall
echo "成功编译generate_data.c文件，生成generate_data可执行程序"
#echo "开始编译elementary_row_transformation.c文件"
#gcc elementray_row_transformation.c -o elementray_row_transformation -g -Wall
#echo "成功编译elementray_row_transformation.c文件，生成elementray_row_transformation可执行程序"
echo "开始编译串行程序lu_decomposition_method.c文件"
gcc lu_decomposition_method.c -o lu_decomposition_method -g -Wall
echo "成功编译串行程序lu_decompostion_method.c文件，生成lu_decomposition_method可执行程序"
echo "开始编译并行程序lu_decompositon_method_parallel.c文件,生成lu_decomposition_method_parallel可执行程序"
mpicc -g -Wall -o lu_decomposition_method_parallel lu_decomposition_method_parallel.c
echo "成功编译并行程序"
echo "开始执行数据生成代码generate_data"
./generate_data
#echo "开始执行elementray_row_transformation"
#./elementray_row_transformation
echo "开始执行lu_decomposition_method"
./lu_decomposition_method
echo "开始执行lu_decomposition_method_parallel"
mpiexec -n 5 ./lu_decomposition_method_parallel
