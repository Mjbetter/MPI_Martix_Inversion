#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>

void file_read(double *num){
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
			num[count++]=tmp;
			tmp=0.0;
			flag=0;
		}
	}
    munmap(addr, sb.st_size);	
    printf("完成文件数据读取\n");
}

int Get_LU(int n,double **num,double **U,double **L){
	int i,j,k;
	/*
	根据公式计算L和U矩阵
	U(k,j) = A(k,j)-∑(t=0,k-1)(L(k,t)*U(t,j))
	L(i,k) = (A(k,j)-∑(t=0,k-1)(L(i,t)*U(t,k)))/U(k,k)
	*/
	/*
	计算U矩阵第一列
	*/
	for(i=0;i<n;++i){
		U[0][i] = num[0][i];
		if(U[0][0]==0)return 0;
	}
	/*
	计算L矩阵的对角线，全为1
	*/
	for(i=0;i<n;++i){
		L[i][i]=1;
	}
	/*
	计算L矩阵第一列
	*/
	for(i=1;i<n;++i){
		L[i][0] = num[i][0]/U[0][0];
	}
	/*
	开始交替计算U和L
	*/
	for(i=1;i<n;++i){
		/*
		计算U
		*/
		double sum1=0;
		for(j=i;j<n;++j){
			sum1=0;
			for(k=0;k<=i-1;++k){
				sum1+=L[i][k]*U[k][j];
			}
			U[i][j]=num[i][j]-sum1;
			if(U[i][i]==0){
//				printf("U[%d][%d]=%.4lf\n",i,i,U[i][i]);
				return i;
			}
		}
		/*
		计算L
		*/
		for(j=i+1;j<=n-1;++j){
			sum1=0;
			for(k=0;k<=i-1;++k){
				sum1+=L[j][k]*U[k][i];
			}
			L[j][i]=(num[j][i]-sum1)/U[i][i];
		}
	}
	return -99999;	
}

int main(){
	struct timeval start,end,start1,end1;
    double seconds;
    double useconds;
    double cost;
    int n=0;
    FILE *fp;
    fp = fopen("data.txt","r");
    fscanf(fp,"%d",&n);
    fclose(fp);
	double *buf=NULL;
	buf = (double*)malloc(sizeof(double)*n*n);
	file_read(buf);
	int i,j,k;
	/*
	动态分配二维数组的内存空间
	LU分解法需要的数组如下：
		1.存储矩阵的二维数组num
		2.存储上三角矩阵的二维数组L
		3.存储下三角矩阵的二维数组U
		4.存储上面三个矩阵的逆矩阵的三个二维数组inverse_num,inverse_L,inverse_U
	*/
	printf("开始分配内存空间\n");
	double **num=(double**)malloc(sizeof(double*)*n);
	double **L=(double**)malloc(sizeof(double*)*n);
	double **U=(double**)malloc(sizeof(double*)*n);
	double **inverse_num=(double**)malloc(sizeof(double*)*n);
	double **inverse_L=(double**)malloc(sizeof(double*)*n);
	double **inverse_U=(double**)malloc(sizeof(double*)*n);
	for(i=0;i<n;++i){
		num[i] = (double*)malloc(sizeof(double)*n);
		L[i] = (double*)malloc(sizeof(double)*n);
		U[i] = (double*)malloc(sizeof(double)*n);
		inverse_num[i] = (double*)malloc(sizeof(double)*n);
		inverse_L[i] = (double*)malloc(sizeof(double)*n);
		inverse_U[i] = (double*)malloc(sizeof(double)*n);
	}
	printf("内存空间分配成功\n");
	/*
	将缓冲区中的数据存入num矩阵
	*/
	printf("开始向num中填入矩阵\n");
	for(i=0;i<n;++i){
		for(j=0;j<n;++j){
			num[i][j] = buf[i*n+j];
		}	
	}
	printf("成功在num中填入矩阵\n");
	/*
	exchange数组用于记录交换记录
	*/
	int *exchange=(int *)malloc(sizeof(int)*n);
	for(i=0;i<n;++i){
		exchange[i]=-99999;
	}
	/*
	从文件中读取矩阵到数组num中
	*/
	
	double sum;
	
	gettimeofday(&start,NULL);
	gettimeofday(&start1,NULL);
	
	/*
	对LU矩阵进行求解,同时检查是否需要进行换行操作
	*/
	int change_row=0;
	while(1){
		change_row=Get_LU(n,num,U,L);
		if(change_row==-99999)break;

		for(i=0;i<n;++i){
			memset(L[i],0,n);
			memset(U[i],0,n);	
		}
		int max_i=change_row;
		double max = num[max_i][change_row];
		/*
		寻找最大的列主元行
		*/
		for(i=change_row+1;i<n;++i){
			if(max<num[i][change_row]&&num[i][change_row]!=0){
				max_i = i;
				max = num[i][change_row];
				break;
			}
		}
		
		if(max_i == change_row){
			break;
		}
		
		printf("第%d行和第%d行发生了换行\n",change_row,max_i);
		double *tmp = num[change_row];
		num[change_row] = num[max_i];
		num[max_i] = tmp;
		exchange[change_row]=max_i;
	}
//	printf("求得的L矩阵为:\n");
//	for(i=0;i<n;++i){
//		for(j=0;j<n;++j){
//			printf("%.4lf ",L[i][j]);
//		}
//		printf("\n");
//	}
//	printf("求的U矩阵为:\n");
//	for(i=0;i<n;++i){
//		for(j=0;j<n;++j){
//			printf("%.4lf ",U[i][j]);
//		}
//		printf("\n");
//	}
	gettimeofday(&end1,NULL);
    seconds  = end1.tv_sec  - start1.tv_sec; // 计算秒数
    useconds = end1.tv_usec - start1.tv_usec; // 计算微秒数
    cost = (seconds * 1000.0) + (useconds / 1000.0); // 计算毫秒数
	printf("计算LU矩阵所耗时间:%.9lf毫秒\n",cost);
	/*
	对U矩阵和L矩阵进行求逆
	存在公式：
	i>j:Linv(i,j)=-Linv(j,j)*∑(i-1,k=j)(L(i,k)*Linv(k,j))
	i<j:Uinv(i,j)=-1/U(i,i)*∑(j,k=i+1)(U(i,k)*Uinv(k,j)))
	*/
	gettimeofday(&start1,NULL);
	for(j=0;j<n;++j){
		for(i=j;i<n;++i){
			if(i==j)inverse_L[i][j]=1/L[i][j];
			else if(i<j) inverse_L[i][j]=0;
			else{
				sum=0;
				for(k=j;k<i;k++){
					sum+=L[i][k]*inverse_L[k][j];
				}
				inverse_L[i][j]=-1*inverse_L[j][j]*sum;
			}
		}
	}
//	printf("求得的L逆矩阵:\n");
//	for(i=0;i<n;++i){
//		for(j=0;j<n;++j){
//			printf("%.4lf ",inverse_L[i][j]);
//		}
//		printf("\n");
//	}
	gettimeofday(&end1,NULL);
    seconds  = end1.tv_sec  - start1.tv_sec; // 计算秒数
    useconds = end1.tv_usec - start1.tv_usec; // 计算微秒数
    cost = (seconds * 1000.0) + (useconds / 1000.0); // 计算毫秒数
	printf("计算L矩阵的逆矩阵所耗时间:%.9lf毫秒\n",cost);
	
	gettimeofday(&start1,NULL);
	for(j=0;j<n;++j){
		for(i=j;i>=0;i--){
			if(i==j)inverse_U[i][j]=1/U[i][j];
			else if(i>j)inverse_U[i][j]=0;
			else{
				sum=0;
				for(k=i+1;k<=j;++k){
					sum+=U[i][k]*inverse_U[k][j];
				}
				inverse_U[i][j]=-1/U[i][i]*sum;
			}
		}
	}
//	printf("求得的U逆矩阵:\n");
//	for(i=0;i<n;++i){
//		for(j=0;j<n;++j){
//			printf("%.4lf ",inverse_U[i][j]);
//		}
//		printf("\n");
//	}
	gettimeofday(&end1,NULL);
    seconds  = end1.tv_sec  - start1.tv_sec; // 计算秒数
    useconds = end1.tv_usec - start1.tv_usec; // 计算微秒数
    cost = (seconds * 1000.0) + (useconds / 1000.0); // 计算毫秒数
	printf("计算U矩阵的逆矩阵所耗时间:%.9lf毫秒\n",cost);
	
	gettimeofday(&start1,NULL);
	for(k=0;k<n;++k){
		for(i=0;i<n;++i){
			double r = inverse_U[i][k];
			for(j=0;j<n;++j){
				inverse_num[i][j]+=r*inverse_L[k][j];
			}
		}
	}
	
	/*
	对换行记录进行遍历，行交换后，逆矩阵的列对应也要交换
	*/
	for(i=0;i<n;++i){
		if(exchange[i]!=-99999){
			printf("第%d列需要和第%d列进行交换\n",i,exchange[i]);
			for(j=0;j<n;++j){
				double tmp = inverse_num[j][i];
				inverse_num[j][i] = inverse_num[j][exchange[i]];
				inverse_num[j][exchange[i]] = tmp;
			}
		}
	}	
	
	gettimeofday(&end1,NULL);
    seconds  = end1.tv_sec  - start1.tv_sec; // 计算秒数
    useconds = end1.tv_usec - start1.tv_usec; // 计算微秒数
    cost = (seconds * 1000.0) + (useconds / 1000.0); // 计算毫秒数
	printf("计算矩阵相乘所耗时间:%.9lf毫秒\n",cost);


	/*
	打印结果代码片段
	*/
	
//	printf("求得的逆矩阵为:\n");
//	for(i=0;i<n;++i){
//		for(j=0;j<n;++j){
//			printf("%.4lf ",inverse_num[i][j]);
//		}
//		printf("\n");
//	}
	
	gettimeofday(&end,NULL);
    seconds  = end.tv_sec  - start.tv_sec; // 计算秒数
    useconds = end.tv_usec - start.tv_usec; // 计算微秒数
    cost = (seconds * 1000.0) + (useconds / 1000.0); // 计算毫秒数
	printf("程序的运行时间为:%.9lf毫秒\n",cost);

	free(exchange);
	for(i=0;i<n;++i){
		free(num[i]);
		free(L[i]);
		free(inverse_L[i]);
		free(U[i]);
		free(inverse_U[i]);
	}		
	free(buf);
	free(num);
	free(L);
	free(U);
	free(inverse_L);
	free(inverse_U);

	return 0;
}

