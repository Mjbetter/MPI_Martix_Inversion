#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>

using namespace std;

int main(){
	int n;
	printf("请输入矩阵的阶数:\n");
	scanf("%d",&n);
	int i,j;
	ofstream f;
	f.open("data.txt",ios::out);
	f << n << endl;
	srand(time(NULL));
	for(i=0;i<n;++i){
		for(j=0;j<n;++j){
			f << (rand()%1000+1)+'0' << " ";
		}
		f << endl;
	}
	f.close();

	return 0;
}
