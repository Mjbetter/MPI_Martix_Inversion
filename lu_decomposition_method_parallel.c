#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*
Exchange函数对矩阵进行换行操作
n表示矩阵维度
m表示需要交换的行
num是二维数据，存储矩阵
buf是num的一维缓冲区
返回被交换的行
*/
int Exchange(int n,int actual_n,int m,double **num,double *buf);

/*
File_Read函数通过文件映射，快速从文件中读取数据，并且对缓冲区中的数据进行格式化处理
num是存储文件数据的缓冲区
*/
void File_Read(double *num);

int main(){
	int rank,size;
	//LU分解法所需要的二维矩阵和阶数n
	double **num=NULL,**L=NULL,**inverse_L=NULL,**inverse_U=NULL;
	double *buf=NULL,*buf1=NULL,*buf_u=NULL, *buf_l=NULL,*buf_inverse_u=NULL,*buf_inverse_l=NULL,*res=NULL,*buffer=NULL;	
	
	/*
	进行MPI的初始化操作
	*/
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	/*
	用于确定程序开始时间和运行结束时间
	*/
	double start_time,finish_time;
	double start,end;
	
	/*
	判断矩阵是否需要扩展，在最后进行矩阵运算的时候，我们需要矩阵维度是进程数量的整数倍
	*/
	int n,actual_n=0;
	
	/*
	进程0需要读取数据，同时将读取到的数据广播到其他进程中	
	*/
 	if(rank==0){
 		printf("这是进程0,我将开始进行数据读取和动态分配存储空间\n");
		FILE *fp;
		fp = fopen("data.txt","r");
		fscanf(fp,"%d",&n);
		fclose(fp);
	}
			
	/*
	广播变量n;
	*/
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	actual_n=n;
	while(actual_n%size!=0)actual_n++;
	/*
	给进程的num数组,buf分配内存空间
	*/
	int i,j,k;
	if(rank==0){
		num=(double**)malloc(sizeof(double*)*actual_n);
		for(i=0;i<actual_n;++i){
			num[i] = (double*)malloc(sizeof(double)*actual_n);
		}		
		buffer = (double*)malloc(sizeof(double)*n*n);		
	}
	buf=(double*)malloc(sizeof(double)*actual_n*actual_n);
	buf1=(double*)malloc(sizeof(double)*actual_n*actual_n);
	if(rank==0){
		int count=0;
		File_Read(buffer);	
		for(i=0;i<actual_n;++i){
			for(j=0;j<actual_n;++j){
				if(i<n&&j<n){
					num[i][j] = buffer[i*n+j];
					buf[count++] = buffer[i*n+j];
				}else{
					num[i][j] = 0;
					buf[count++] = 0;
				}
			}
		}
	}
	
	/*
	开辟进程所需要的进程空间
	*/
	
	buf_inverse_u=(double*)malloc(sizeof(double)*actual_n*actual_n);	
	buf_inverse_l=(double*)malloc(sizeof(double)*actual_n*actual_n);
	buf_u=(double*)malloc(sizeof(double)*actual_n*actual_n);	
	buf_l=(double*)malloc(sizeof(double)*actual_n*actual_n);

	if(rank==0)L=(double **)malloc(sizeof(double*)*actual_n);
	inverse_U=(double **)malloc(sizeof(double*)*actual_n);
	inverse_L=(double **)malloc(sizeof(double*)*actual_n);	
	for(i=0;i<actual_n;++i){
		if(rank==0)L[i] = (double*)malloc(sizeof(double)*actual_n);
		inverse_U[i] = (double*)malloc(sizeof(double)*actual_n);
		inverse_L[i] = (double*)malloc(sizeof(double)*actual_n);
	}
	/*
	开辟用于记录交换记录的数组
	*/
	int *exchange_array=NULL;
	exchange_array = (int*)malloc(sizeof(int)*actual_n);
	for(i=0;i<actual_n;++i){
		exchange_array[i] = -99999;
	}
	if(rank==0){
		printf("完成内存分配\n");
	}
	/*
	计时开始
	*/
	start_time = MPI_Wtime();
	
	/*
	L矩阵和U矩阵的计算相互依赖，如果使用并行来计算LU矩阵，那么会让进程间的通信消耗巨大，所以我们直接采用串行的方式计算出LU矩阵
	对LU矩阵进行求解,同时检查是否需要进行换行操作
	*/
	start = MPI_Wtime();

	/*
	现在每个进程都有矩阵num的数据，L矩阵的数据，U矩阵的数据
	我们可以将L矩阵的逆，U矩阵的逆通过大量的进程同时进行计算
	*/
	int rows_pre_process = actual_n / size;
	int start_row = rank*rows_pre_process;
	double *max_line = (double*)malloc(sizeof(double)*actual_n);
	double *local_num = (double*)malloc(sizeof(double)*rows_pre_process*actual_n);
	
	/*
	用列主元消元法求解U矩阵，同时可以计算出L矩阵
	*/
	start = MPI_Wtime();
	while(1){
		int flag=0;
		MPI_Scatter(buf,rows_pre_process*actual_n,MPI_DOUBLE,local_num,rows_pre_process*actual_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
		for(i=0;i<actual_n;++i){	
			if(rank==0){
				if(buf[i*actual_n+i]==0&&i<n){
					/*
					恢复buf缓冲数组
					*/
					flag=1;
					for(j=0;j<actual_n;++j){
						for(k=0;k<actual_n;++k){
							buf[j*actual_n+k] = num[j][k];
						}
					}
					exchange_array[i]=Exchange(n,actual_n,i,num,buf);
				}
				if(flag==0){
					for(j=0;j<actual_n;++j){
						max_line[j] = buf[i*actual_n+j];
					}
					/*
					每次进行消元前，可以计算出L矩阵的第i列
					*/
					for(j=i;j<actual_n;++j){
						if(i==j)L[j][i] = 1;
						else L[j][i] = buf[j*actual_n+i] / buf[i*actual_n+i];
					}	
				}	
			}
			MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			if(flag==1)break;
			MPI_Bcast(max_line,actual_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
			/*
			每个进程会分配行数均衡的数据，然后我们将这些行和主元行进行消元，然后聚集回主进程，等待下一次分配
			*/
			for(j=0;j<rows_pre_process;++j){
				if(i<j+start_row){
					double tmp=max_line[i];
					double mutiple_num = local_num[j*actual_n+i]/tmp;
					if(local_num[j*actual_n+i]+mutiple_num*tmp==0){
						for(k=i;k<actual_n;++k){
							local_num[j*actual_n+k] += mutiple_num*max_line[k];
						}
					}
					else{
						for(k=i;k<actual_n;++k){
							local_num[j*actual_n+k] -= mutiple_num*max_line[k];
						}					
					}
						
				}
			}
			/*
			结果聚集到num数组上，最后num数据就会变成U矩阵
			*/
			MPI_Gather(local_num,rows_pre_process*actual_n,MPI_DOUBLE,buf,rows_pre_process*actual_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
		}
		if(flag==0)break;		
	}
	if(rank==0){
		for(i=0;i<actual_n;++i){
			for(j=0;j<actual_n;++j){
					buf_u[i*actual_n+j] = buf[i*actual_n+j];
					buf_l[i*actual_n+j] = L[i][j];
				
			}
		}
	}
	if(rank==0){
		end = MPI_Wtime();
		printf("LU计算消耗时间%.9lf\n",(end-start)*1000.0);
	}
	start = MPI_Wtime();
	MPI_Bcast(buf_u, actual_n*actual_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(buf_l, actual_n*actual_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	

	if(rank==0){
		end = MPI_Wtime();
		printf("LU广播通信消耗时间%.9lf\n",(end-start)*1000.0);
	}	
	
	/*
	local_x作为存储本进程计算结果 
	*/
	double *local_u = (double*)malloc(sizeof(double)*actual_n*rows_pre_process);
	double *local_l = (double*)malloc(sizeof(double)*actual_n*rows_pre_process);
	
	/*
	对U矩阵进行求逆
	*/
	start = MPI_Wtime();
	int count=0;
	for(j=start_row;j<start_row+rows_pre_process;++j){
		for(i=j;i>=0;i--){
			if(i==j)inverse_U[i][j]=1/buf_u[i*actual_n+j];
			else if(i>j)inverse_U[i][j]=0;
			else{
				double sum=0;
				for(k=i+1;k<=j;++k){
					sum+=buf_u[i*actual_n+k]*inverse_U[k][j];
				}
				inverse_U[i][j]=-1/buf_u[i*actual_n+i]*sum;
			}
		}
	}
	for(j=start_row;j<start_row+rows_pre_process;++j){
		for(i=0;i<actual_n;++i){
			local_u[count++] = inverse_U[i][j];
		}
	}

	if(rank==0){
		end = MPI_Wtime();
		printf("计算U逆矩阵耗时%.9lf毫秒\n",(end-start)*1000.0);
	}
	/*
	将各个进程的计算结果聚集到0号进程的buf_inverse_u数组中
	*/
	if(rank==0){
		printf("开始U逆矩阵的数据聚合\n");
	}
	start = MPI_Wtime();
	
	MPI_Gather(local_u,rows_pre_process*actual_n,MPI_DOUBLE,buf,rows_pre_process*actual_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if(rank==0){
		end = MPI_Wtime();
		printf("聚合U逆矩阵耗时%.9lf毫秒\n",(end-start)*1000.0);
	}	
	count=0;
	start = MPI_Wtime();
	for(j=start_row;j<start_row+rows_pre_process;++j){
		for(i=j;i<actual_n;++i){
			if(i==j)inverse_L[i][j]=1/buf_l[i*actual_n+j];
			else if(i<j) inverse_L[i][j]=0;
			else{
				double sum=0;
				for(k=j;k<i;k++){
					sum+=buf_l[i*actual_n+k]*inverse_L[k][j];
				}
				inverse_L[i][j]=-1*inverse_L[j][j]*sum;
			}
		}
	}
	for(j=start_row;j<start_row+rows_pre_process;++j){
		for(i=0;i<actual_n;++i){
			local_l[count++] = inverse_L[i][j];
		}
	}
	if(rank==0){
		end = MPI_Wtime();
		printf("计算L逆矩阵耗时%.9lf毫秒\n",(end-start)*1000.0);
	}	
	//将各个进程的计算结果聚集到0号进程的buf_inverse_l数组中
	start = MPI_Wtime();
	if(rank==0){
		printf("开始L逆矩阵的数据聚合\n");
	}
	MPI_Gather(local_l,rows_pre_process*actual_n,MPI_DOUBLE,buf1,rows_pre_process*actual_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if(rank==0){
		end = MPI_Wtime();
		printf("聚合L逆矩阵耗时%.9lf毫秒\n",(end-start)*1000.0);
	}	
	/*
	将矩阵顺序进行整理
	*/
	for(i=0;i<n;++i){
		for(j=0;j<n;++j){
			buf_inverse_u[i*actual_n+j] = buf[j*actual_n+i];
			buf_inverse_l[i*actual_n+j] = buf1[j*actual_n+i];
		}
	}	
	
	MPI_Bcast(buf_inverse_l, actual_n*actual_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	/*
	每个进程计算一部分的数据，然后通过聚集将集中到进程0
	res用于存储最终结果，一维数组便于进行进程间通信
	*/
	
	/*
	为计算所需要的内存空间进行分配
	将3号进程中的U逆矩阵存储在一维数组中进行分发，然后多个进程同时进行矩阵的乘法运算加快速度
	*/
	start = MPI_Wtime();
	double *local_inverse_u = (double*)malloc(sizeof(double)*rows_pre_process*actual_n);
	start = MPI_Wtime();
	MPI_Scatter(buf_inverse_u,rows_pre_process*actual_n,MPI_DOUBLE,local_inverse_u,rows_pre_process*actual_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if(rank==0){
		end = MPI_Wtime();
		printf("分发U逆矩阵数据给各个进程消耗%.9lf毫秒\n",(end-start)*1000.0);
	}	
	double *local_res = (double*)malloc(sizeof(double)*rows_pre_process*actual_n);
	for(k=0;k<actual_n;++k){
		for(i=0;i<rows_pre_process;++i){
			double r = local_inverse_u[i*actual_n+k];
			for(j=0;j<actual_n;++j){
				local_res[i*actual_n+j] += r*buf_inverse_l[k*actual_n+j];
			}
		}
	}
	if(rank==0){
		end = MPI_Wtime();
		printf("L逆矩阵和U逆矩阵相乘耗时%.9lf毫秒\n",(end-start)*1000.0);
	}
	/*
	将所有进程的计算结果全部聚集到进程0中，并进行打印
	*/
	if(rank==0){
		res=(double*)malloc(sizeof(double)*actual_n*actual_n);	
		printf("即将聚集结果数据\n");
	}
	start = MPI_Wtime();
	MPI_Gather(local_res,rows_pre_process*actual_n,MPI_DOUBLE,res,rows_pre_process*actual_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if(rank==0){
		end = MPI_Wtime();
		printf("聚合相乘结果耗时%.9lf毫秒\n",(end-start)*1000.0);
	}
	if(rank==0){
		for(i=0;i<n;++i){
			if(exchange_array[i]!=-99999){
				printf("第%d列需要和第%d列进行交换\n",i,exchange_array[i]);
				for(j=0;j<n;++j){
					double tmp = res[j*actual_n+i];
					res[j*actual_n+i] = res[j*actual_n+exchange_array[i]];
					res[j*actual_n+exchange_array[i]] = tmp;
				}
			}
		}		
	}
	
	/*
	打印结果代码片段
	*/
//	if(rank==0){
//		printf("程序运行结束,求得的逆矩阵为:\n");
//		int i,j;
//		for(i=0;i<n;++i){
//			for(j=0;j<n;++j){
//				printf("%.4lf ",res[i*actual_n+j]);
//			}
//			printf("\n");
//		}
//	}
	
	finish_time=MPI_Wtime();
	if(rank==0){
		printf("程序所需时间为:%.9lf毫秒\n",(finish_time-start_time)*1000.0);
	}
	/*
	进程资源销毁
	*/		
	free(max_line);
	free(local_num);
	if(rank==0)free(buffer);
	free(buf);	
	free(buf1);
	free(buf_inverse_u);
	free(buf_inverse_l);
	free(buf_u);
	free(buf_l);
	free(local_inverse_u);
	free(local_res);
	free(local_l);
	free(local_u);
	if(rank==0){
		free(res);
	}
	free(exchange_array);
	for(i=0;i<actual_n;++i){
		if(rank==0)free(L[i]);
		if(rank==0)free(num[i]);
		free(inverse_L[i]);
		free(inverse_U[i]);
	}		
	if(rank==0)free(num);
	if(rank==0)free(L);
	free(inverse_L);
	free(inverse_U);
	
	
	
	MPI_Finalize();
	return 0;
}


int Exchange(int n,int actual_n,int m,double **num,double *buf){
	int i;
	double max=buf[m*n+m];
	int max_i = m;
	for(i=m+1;i<n;++i){
		if(max<buf[i*actual_n+m]&&buf[i*actual_n+m]!=0){
			max = buf[i*actual_n+m];
			max_i = i;
		}
	}
	int j;
	for(i=0;i<actual_n;++i){
		double tmp = buf[m*actual_n+i];
		buf[m*actual_n+i] = buf[max_i*n+i];
		buf[max_i*n+i] = tmp;
	}
	for(i=0;i<n;++i){
		for(j=0;j<n;++j){
			num[i][j] = buf[i*n+j];
		}
	}	
	return max_i;
}

void File_Read(double *buffer){
	printf("开始从文件中读取数据\n");
    int fd;
    char *addr;
    struct stat sb;
	
    fd = open("data.txt", O_RDWR);
    fstat(fd, &sb);
    /*
	将fd句柄指向的文件映射到addr内存地址
	*/
    addr = mmap(NULL, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);
	
	/*
	对文件内的数据进行格式化处理
	*/
	long long size = strlen(addr);
	long long i,start=0;
	while(addr[start]!='\n')start++;
	long long count=0;
	long long flag=0;
	double tmp=0.0;
	for(i=start;i<size;++i){
		if(addr[i]=='\n')i++;
		if(addr[i]>='0'&&addr[i]<='9'){
			if(flag==0)tmp=tmp*10+(addr[i]-'0');
			else{
				flag*=10;
				tmp=tmp+(addr[i]-'0')/flag;
			}
		}
		else if(addr[i]=='.')flag=1;
		else if(addr[i]==' '&&addr[i+1]!=' '){
			buffer[count++]=tmp;
			tmp=0.0;
			flag=0;
		}
	}
    munmap(addr, sb.st_size);	
    printf("完成文件数据读取\n");
}
