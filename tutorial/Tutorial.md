# 人工智能项目入手向导

本说明主要介绍一些入手人工智能相关项目所需的必要`软件`及`知识`储备，用于帮助新人快速入门或上手相关项目。如果你想更加全面的了解人工智能项目开发相关的技术，也可以根据本说明提到的一些关键词自行搜索学习。

## 目录

- [软件安装](#软件安装)
- [Python 环境配置](#python-环境配置)
- [项目开发相关知识](#项目开发相关知识)
  - [Python 编程语言](#python-编程语言)
  - [Git 相关](#git-相关)
  - [镜像资源修改](#镜像资源修改)

## 软件安装

### 集成开发环境（IDE）配置

鉴于现阶段人工智能项目主要基于[Python](https://www.python.org/)编程语言进行开发，我们推荐入门人员首选安装[Pycharm](https://www.jetbrains.com/pycharm/)并配置[Conda](https://www.anaconda.com/products/distribution)环境使用。另外，除了`Python`以外，项目开发过程中还会使用其他的编程语言（如`C++`，`Javascript`等）、标记语言（如`XML`、`HTML`等）以及特殊的文件格式（如`Json`等），而[Vscode](https://code.visualstudio.com/)因为有丰富的插件库用以支持上述语言及文件格式，所以我们也推荐新手安装[Vscode](https://code.visualstudio.com/)。

下面将就`Pycharm`，`Conda`以及`VScode`的安装及使用进行简单的说明。

#### Pycharm 的安装

浏览器输入 https://www.jetbrains.com/pycharm/ 进入 Pycharm 的下载页面，根据电脑系统选择下载 Windows，MACOS 还是 Linux 版本。

![img.png](img.png)

**Note**

> 安装 Pycharm 社区版即可

#### Anaconda 的安装

浏览器输入 https://www.anaconda.com/products/distribution 进入 Anaconda 的下载页面，点击`Download`会自动下载电脑适配版本。

![img_1.png](img_1.png)

**Important**

> 在安装 Anaconda 过程中，注意勾选将 Anaconda 加入 PATH 环境变量，安装结束之后可查看系统环境变量确认

<div align=center><img src="img_14.png"></div>

> 系统环境变量查找路径：此电脑-属性（右键单击）-高级系统设置-环境变量-Path 环境变量（双击可查看）

<div align=center><img width=500 src="img_3.png"/></div>

#### VScode 安装 (可选，建议安装)

浏览器输入 https://code.visualstudio.com/ 进入 VScode 的下载页面，单击`Download for Windows`即可下载 Windows 版本的 VScode。

![img_2.png](img_2.png)

### 爬虫相关组件

人工智能项目开发中涉及到利用爬虫开展数据挖掘的任务时，会使用[Selenium](https://www.selenium.dev/)等框架控制浏览器开展数据爬取，因此需要安装浏览器及对应的自动控制程序。

本说明主要介绍`Chrome`浏览器及对应的`ChromeDriver`组件的安装，其他的浏览器（如`Edge`、`Safari`）也可以按类似的方法进行进行安装配置。

#### Chrome 浏览器

浏览器输入 https://www.google.cn/chrome/index.html 进入 Chrome 浏览器下载页面，点击`下载Chrome`即可下载系统适配的 Chrome 浏览器。

![img_4.png](img_4.png)

#### ChromeDriver 组件

ChromeDriver 组件需要和对应版本的 Chrome 浏览器一起搭配使用，因此安装完`Chrome`浏览器之后，点击`浏览器-设置-关于Chrome`查看安装的 Chrome 版本，例如我当前的 Chrome 版本号为`104.0.5112.82`。

![img_5.png](img_5.png)

确认安装的`Chrome`版本后，浏览器输入 http://chromedriver.storage.googleapis.com/index.html 进入`ChromeDriver`下载页面，根据 Chrome 版本号选择 ChromeDriver 版本，如果没有对应的版本，选择最相近的版本，例如我这里选择`104.0.5112.79`文件夹进入。

![img_6.png](img_6.png)

进入`104.0.5112.79`之后，即可选择对应系统的`ChromeDriver`进行下载安装。

![img_7.png](img_7.png)

**note**

> 所谓的 ChromeDriver 其实就是一个压缩包，里面存放有 chromedriver.exe（Win32 版本）文件。

**important**

> 为了在 Python 环境下使用 ChromeDriver，请将 chromedriver.exe 的位置加入到上述提到的 PATH 环境变量中（另一个更为简单的方法是把该 exe 文件直接放到 conda 目录下的 Scripts 文件夹中（已加入 PATH 环境变量），例如我计算机的位置：D:\Anaconda\Scripts）。

## Python 环境配置

下面将就如何利用 conda 配置 Python 环境进行说明，首先介绍 PowerShell 中 conda 的配置，然后讲述 conda 的一些常用命令，最后介绍如何利用 Pycharm 配置 Conda 环境。

### PowerShell 中配置 Conda（Windows 平台）

对于 Windows 平台，初次使用 conda 时，首先在任务栏的搜索框中搜索`PowerShell`，点击`以管理员身份运行`

![img_8.png](img_8.png)

然后输入以下命令（目的是将 powershell 和 anaconda 进行关联）：

```
Set-ExecutionPolicy RemoteSigned
```

会显示如下界面：

![img_9.png](img_9.png)

选择`Y`进行确定，之后就可以在 PowerShell 中正常使用 conda 了。

### Conda 创建新环境

启动 conda 时，默认显示的是 base 环境，即上图第二行最左侧的`(base)`。

对于项目开发来说，我们推荐每一个项目都在一个新的环境下开发，这可以极大的避免项目之间因为所需的 Python 以及 Python 包的版本不同而产生冲突。

conda 中创建新环境（例如 python38）可以使用下述命令：

```
conda create -n python38
```

创建好`python38`环境后，如果想查看当前 conda 的所有环境，可以使用下述命令：

```
conda env list
```

![img_10.png](img_10.png)

可以看到，`python38`环境已经创建成功啦。

### Conda 激活环境

当我们想使用新的`python38`环境时，我们需要先激活该环境，使用下述命令：

```
conda activate python38
```

激活后，`PS`符号左侧的`(base)`就会变成`(python38)`。

### 使用 Conda 安装 Python 及各种包

对于新创建的环境默认是没有 Python 的，所以我们首先需要安装 Python，利用 Conda 安装 Python 时我们可以同时指定 Python 的版本，使用下述命令：

```
conda install python=3.8
```

使用下述命令查看 python 是否安装成功：

```
conda list
```

![img_11.png](img_11.png)

可以看到，版本号为`3.8.13`的 python 已经成功安装在新环境下；另外，除了 python 之外，conda 还帮助我们安装了其他安装 python 所需的依赖（如`setuptools`等）。

此外，`conda install`命令不仅可以用来安装 python，其他各种 python 包也可以使用该命令进行安装，只需将 python 改为包的名字即可。

### Pycharm 中使用 Conda

在日常的项目开发中，我们更多的是直接利用 Pycharm 启动一个 conda 环境进行开发，需要做的就是在首次打开项目时进行简单的配置：

例如，当我们利用打开[GVasp](https://github.com/Rasic2/gvasp)项目时， 点击右下角的`Python解释器`，如图中右下角红色框的位置：

![img_12.png](img_12.png)

点击之后，选择`添加解释器-Conda环境-现有环境`，如图：

![img_13.png](img_13.png)

在`解释器`位置可以更改你想修改的 conda 环境，选择之后点击确认，项目的 conda 环境就配置好了。

配置好 conda 环境后，如果后续想进行 conda 包的安装或环境的更改，可以直接利用 Pycharm 自带的 PowerShell 窗口进行操作，不用在搜索 PowerShell 启动使用了。

**note**

> 在 Pycharm 中使用 Powershell 可以点击 Pycharm 底部菜单栏的终端启动 Powershell（小箭头选择）

## 项目开发相关知识

### Python 编程语言

### Git 相关

下面我们将简单介绍 git 相关的命令以方便新手入门如何使用 github 进行项目代码的管理。

对于项目开发来说，往往涉及到多人合作开发一个项目以及项目代码开发过程中的版本更新迭代，于是需要一种工具来管理代码开发的整个过程，涉及代码的修改历史，版本回溯，分支合并等等，于是 git 就产生了，而 github 则是一个代码托管的平台，更多有关 git 及 github 的知识感兴趣的可以自行搜索学习。

常用的 git 命令如下：

- 下载别人的项目代码

```
git clone https://github.com/Rasic2/gvasp.git ./gvasp
```

其中 `https:*` 代表你想下载的项目代码，`./gvasp` 代表你想把代码下载到本地的哪个位置。

- 初始化自己的仓库

```
git init
```

- 查看当前项目仓库的配置

```
git config -l
```

- 将当前目录下的所有文件加入历史追踪

```
git add .
```

- 将当前所作的修改保存为一次本地版本更新

```
git commit -m "comment"
```

其中 `comment` 表示你对本次修改的一个批注。

- 查看当前项目代码的历史日志

```
git log
```

- 给远程仓库起一个别名

```
git remote add origin git@github.com:Rasic2/gvasp.git
```

其中 `origin`代表远程仓库地址的一个别名，`git@github*`代表远程仓库的地址。

- 将本地的最新版本同步到远程仓库

```
git push origin master
```

其中 `origin`代表远程仓库地址的一个别名，`master`代表一个代码分支的名字。

- 将远程仓库的最新代码同步到本地

```
git pull origin master
```

由于 git 相关的知识较为丰富，本说明仅列了较为常用的 git 命令，更多相关的知识可以参考下述链接学习。

#### 更多参考资源

[1. 常用 Git 命令总结](https://zhuanlan.zhihu.com/p/384819351)

[2. Git 教程](https://www.runoob.com/git/git-tutorial.html)

[3. Git reset 命令的使用](https://www.jianshu.com/p/cbd5cd504f14)

[4. 细读 Git | 让你弄懂 origin、HEAD、FETCH_HEAD 相关内容](https://developer.aliyun.com/article/919354)

[5. github fork 与pull request](https://blog.csdn.net/benben0729/article/details/83031135)

### 镜像资源修改

#### Conda 镜像修改

对于 conda，由于国外资源连通率太低，我们通常会使用 conda 的国内镜像资源，方法是在`Users`目录下（如 C:\Users\hui_zhou）新建一个`.condarc`文件，并修改如下（使用清华镜像）：

```
show_channel_urls: true
default_channels:
  - http://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - http://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - http://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch/
```

#### Pip 镜像修改

除了 conda 之外，对于少数包，如果使用 pip 进行安装时，也许使用国内的镜像资源进行加速，方法是在用户目录下新建`.pip/pip.conf`文件（先建一个`.pip`文件夹），然后修改如下（使用阿里镜像）：

```
[global]
index-url = https://mirrors.aliyun.com/pypi/simple
```

## 贡献者

[Hui Zhou](https://github.com/Rasic2), [DongZhi Li](https://github.com/mastreina)

## 版权

本说明的版权归[GongGroup](https://github.com/GongGroup)所有，更新维护工作暂由[Hui Zhou](https://github.com/Rasic2)及[DongZhi Li](https://github.com/mastreina)负责。
