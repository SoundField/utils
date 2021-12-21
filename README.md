[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)
[![license](https://img.shields.io/npm/l/debug)](https://github.com/SoundField/utils/LICENSE)
# 一些实用绘图代码和程序软件

目录和内容格式：  
```
>[坚果云](##坚果云)  
## 坚果云
**软件（代码）用途：**  
多平台同步软件，效率高，开源免费。  
**使用方法：**    
官网：https://www.jianguoyun.com/s/downloads  
yay -S nutstore  
**注意事项：**  
yay -S python-gobject#如果出现坚果云打不开或者没有登录页面的情况，安装坚果云相关依赖
Linux出现客户端白屏已经解决了，但是需要等AUR更新。
```
# 目录
[Matlab绘图](#Matlab)  
>[声场伪彩图](##伪彩图)  

[Windows软件](#Windows)  
>[TeamViews](##TeamsViews)  

[Latex软件](#Latex)  
>[坚果云](##坚果云)  

# Matlab

## 伪彩图
**软件（代码）用途：**  
绘制声压传播损失伪彩图，有两个版本  
具体 [说明](./Codes/Matlab/Rams伪彩图/MATLAB读取程序/说明.txt)  
**使用方法：**    
文件目录在这里 [here](./Codes/Matlab/Rams伪彩图)   
Rams两层海底文件夹中有两个版本的可执行程序，附上了绘图程序在Matlab读取程序中；   
首先执行可执行程序，计算得到的结果，纯声压版本会快很多。  
然后Matlab改写添加计算结果的路径并执行
```(matlab)
  readdata_pcolor.m （均匀网格）
  # readdata_pcolor2.m （非均匀网格）
```  
**注意事项：**   
需要改写路径  

# Windows

## TeamViews  
**软件（代码）用途：**   
远程访问、远程支持、跨平台、跨设备，免费  
**使用方法：**     
https://www.teamviewer.cn/cn/  
**注意事项：**   
要注册账号，选远程控制（使用密码）可以直接访问，另一端需要打开TeamViews。





# Latex

## 坚果云
**软件（代码）用途：**  
多平台同步软件，效率高，开源免费。  
**使用方法：**    
官网：https://www.jianguoyun.com/s/downloads  
yay -S nutstore  
**注意事项：**  
yay -S python-gobject#如果出现坚果云打不开或者没有登录页面的情况，安装坚果云相关依赖
Linux出现客户端白屏已经解决了，但是需要等AUR更新。



## License

[MIT © SoundField.](../LICENSE)
