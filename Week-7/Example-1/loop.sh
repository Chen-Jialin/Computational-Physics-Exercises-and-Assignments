#! /bin/bash
# 本文件用于循环改变晶格常数进行计算以获得能量最低情况下的晶格常数
BIN=$HOME/bin/vasp_std # 设置VASP软件的路径
rm WAVECAR SUMMARY.fcc # 删除上次计算生成的波函数文件（因为是从晶体结构开始算起）和上次运行本脚本时生成的日志文件
for i in  3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 ; do # 将所列的这串数字一次代入i并做以下的事情
cat >POSCAR <<!
fcc:
   $i
 0.5 0.5 0.0
 0.0 0.5 0.5
 0.5 0.0 0.5
   1
cartesian
0 0 0
!# 将下列内容写入POSCAR文件，到!处结束
echo "a= $i" ; $BIN # 向屏幕输出本次循环所用的晶格常数并调用VASP软件开始运算
E=`awk '/F=/ {print $0}' OSZICAR` ; echo $i $E  >>SUMMARY.fcc # 将OSZICAR中'F='这个字符串及其之后的内容赋给E这个变量，将i和E这两个变量输入日志文件
done # 本次循环做完了
cat SUMMARY.fcc # 向屏幕输出本次脚本运行生成的日志文件
