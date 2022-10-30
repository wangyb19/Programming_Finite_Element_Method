# Programming the Finite Element Method

This repo contains the source code of *Programming the Finite Element Method Fifth Edition*, I. M. Smith, D. V. Griffiths and L. Margetts.

## 编译-链接

以编译链接运行本书第一个程序*p41.f03*为例

1. 将[远程仓库](https://github.com/wangyb19/Programming_Finite_Element_Method)下载到本地：`git clone https://github.com/wangyb19/Programming_Finite_Element_Method`

2. 进入项目library文件夹 `cd ./Programming_Finite_Element_Method/library/ `
3. 进入main文件夹 `cd main`/
4. 编译所有子程序和函数以及模块文件：`gfortran -c *.f03`

5. 使用 *ar* 创建静态库: `ar -r mainlib.a *.o`
6. 删除中间文件： `rm *.o`
7. 将静态库文件和模块文件复制到和*p41.f03*同一目录下
   * `cp ./mainlib.a ../../chap04/`
   * `cp ./main.mod ../../chap04/`
8. 进入geom文件夹：`cd ../geom/`
9. 编译所有子程序和函数以及模块文件：`gfortran -c *.f03`
10. 使用 *ar* 创建静态库: `ar -r geomlib.a *.o`
11. 删除中间文件： `rm *.o`
12. 将静态库文件和模块文件复制到和*p41.f03*同一目录下
    * `cp ./geomlib.a ../../chap04/`
    * `cp ./geom.mod ../../chap04/`
13. 将另一个静态库文件*arpacklib.a*也复制到和*p41.f03*同一目录下
    * `cp ../arpacklib.a ../../chap04/`

5. 进入*p41.f03*所在文件夹: `cd ../../chap04`
6. 对*p41.f03*编译链接生成可执行程序*p41*
   * `gfortran -o p41 p41.f03 arpacklib.a geomlib.a mainlib.a`
7. 运行可执行程序并按照提示输入相关信息进行计算: `./p41`