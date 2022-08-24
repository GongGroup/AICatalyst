# 人工智能项目教程

本教程主要介绍一些入手人工智能相关项目所需的必要`软件`及`知识`储备，用于帮助新人快速入门或上手相关项目。如果你想更加全面的了解人工智能项目开发相关的技术，也可以根据本教程提到的一些关键词自行搜索学习。

## 目录

- [软件安装](#软件安装)
- [Python 环境配置](#python-环境配置)
- [项目开发相关知识](#项目开发相关知识)
  - [Git 相关](#git-相关)
  - [镜像资源修改](#镜像资源修改)
- [学习资料推荐](#学习资料推荐)
  - [Python 编程语言](#Python-编程语言)
  - [Python 数据分析](#Python-数据分析)
  - [Python 网络爬虫](#Python-网络爬虫)
  - [Python 机器学习](#Python-机器学习)
  - [拓展书目](#拓展书目)
- [常用库的在线手册](#常用库的在线手册)

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

浏览器输入 https://www.anaconda.com/products/distribution 进入 Anaconda 的下载页面，点击`Download`会自行下载电脑适配版本。

![img_1.png](img_1.png)

**Important**

> 在安装 Anaconda 过程中，注意勾选将 Anaconda 加入 PATH 环境变量，安装结束之后可查看系统环境变量确认

<div align=center><img width=400 src="img_14.png"></div>

> 系统环境变量查找路径：此电脑-属性（右键单击）-高级系统设置-环境变量-Path 环境变量（双击可查看）

<div align=center><img width=400 src="img_3.png"/></div>

#### VScode 安装 (可选，建议安装)

浏览器输入 https://code.visualstudio.com/ 进入 VScode 的下载页面，单击`Download for Windows`即可下载 Windows 版本的 VScode。

![img_2.png](img_2.png)

### VSCode 插件安装

与 Pycharm 不同，VScode 的强大主要得益于其具有丰富且优秀的插件来帮助开发者更便利的开展项目开发。

在 VSCode 中进行插件安装特别简单，只需要在`插件搜索栏`中输入你想安装的插件名字，即可进行安装。`插件搜索栏`如图中红框所示：

![img_15.png](img_15.png)

推荐安装的 VSCode 插件：

- [C/C++: 官方 C/C++语言插件](https://github.com/Microsoft/vscode-cpptools.git)
- [Python: 官方 Python 语言插件](https://github.com/Microsoft/vscode-python)
- [Prettier: 代码格式化插件](https://prettier.io/)
- [Markdown Preview Enhanced: Markdown 预览插件](https://shd101wyy.github.io/markdown-preview-enhanced/#/)
- [Git Graph: 查看 Git 提交历史插件](https://marketplace.visualstudio.com/items?itemName=mhutchie.git-graph)
- [Live Server: 本地服务器插件](https://github.com/ritwickdey/vscode-live-server.git)
- [Code Spell Checker: 检查拼写插件](https://github.com/streetsidesoftware/vscode-spell-checker.git)

### 爬虫相关组件

人工智能项目开发中涉及到利用爬虫开展数据挖掘的任务时，会使用[Selenium](https://www.selenium.dev/)等框架控制浏览器开展数据爬取，因此需要安装浏览器及对应的自动控制程序。

本说明主要介绍`Chrome`浏览器及对应的`ChromeDriver`组件的安装，其他的浏览器（如`Edge`、`Safari`等）也可以按类似的方法进行进行安装配置。

#### Chrome 浏览器

浏览器输入 https://www.google.cn/chrome/index.html 进入 Chrome 浏览器下载页面，点击`下载Chrome`即可下载系统适配的 Chrome 浏览器。

![img_4.png](img_4.png)

#### ChromeDriver 组件

ChromeDriver 组件需要和对应版本的 Chrome 浏览器一起搭配使用，因此安装完`Chrome`浏览器之后，点击`浏览器-设置-关于Chrome`查看安装的 Chrome 版本，例如我当前的 Chrome 版本号为`104.0.5112.82`。

![img_5.png](img_5.png)

确认安装的`Chrome`版本后，浏览器输入 http://chromedriver.storage.googleapis.com/index.html 进入`ChromeDriver`下载页面，根据 Chrome 版本号选择 ChromeDriver 版本，如果没有对应的版本，选择最相近的版本，例如我这里选择`104.0.5112.79`文件夹进入。

![img_6.png](img_6.png)

进入`104.0.5112.79`之后，即可选择对应系统的`ChromeDriver`进行下载安装。

<div align=center><img width=600 src="img_7.png"/></div>

**note**

> 所谓的 ChromeDriver 其实就是一个压缩包，里面存放有 chromedriver.exe（Win32 版本）文件。

**important**

> 为了在 Python 环境下使用 ChromeDriver，请将 chromedriver.exe 的位置加入到上述提到的 PATH 环境变量中（另一个更为简单的方法是把该 exe 文件直接放到 conda 目录下的 Scripts 文件夹中（已加入 PATH 环境变量），例如我计算机的位置：D:\Anaconda\Scripts）。

## Python 环境配置

下面将就如何利用 conda 配置 Python 环境进行说明，首先介绍 PowerShell 中 conda 的配置，然后讲述 conda 的一些常用命令，最后介绍如何利用 Pycharm 配置 Conda 环境。

### PowerShell 中配置 Conda（Windows 平台）

对于 Windows 平台，初次使用 conda 时，首先在任务栏的搜索框中搜索`PowerShell`，点击`以管理员身份运行`

<div align=center><img width=500 src="img_8.png"/></div>

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

<div align=center><img width=400 src="img_10.png"/></div>

可以看到，`python38`环境已经创建成功啦。

### Conda 激活环境

当我们想使用新的`python38`环境时，需要先激活该环境，使用下述命令：

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

<div align=center><img width=500 src="img_11.png"/></div>

可以看到，版本号为`3.8.13`的 python 已经成功安装在新环境下；另外，除了 python 之外，conda 还帮助我们安装了其他安装 python 所需的依赖（如`setuptools`等）。

此外，`conda install`命令不仅可以用来安装 python，其他各种 python 包也可以使用该命令进行安装，只需将 `python` 改为包的名字即可（如`numpy`, `PyQuery`等）。

### Pycharm 中使用 Conda

在日常的项目开发中，我们更多的是直接利用 Pycharm 启动一个 conda 环境进行开发，需要做的仅是在首次打开项目时进行简单的配置：

以[GVasp](https://github.com/Rasic2/gvasp)项目为例，当我们在 Pycharm 中打开该项目时， 点击右下角的`Python解释器`，如图中右下角红色框的位置：

![img_12.png](img_12.png)

点击之后，选择`添加解释器-Conda环境-现有环境`，如图：

![img_13.png](img_13.png)

在`解释器`位置更改你想修改的 conda 环境，选择之后点击确认，项目的 conda 环境就配置好了。

配置好 conda 环境后，如果后续想进行 conda 包的安装，可以直接利用 Pycharm 自带的 PowerShell 窗口进行操作，无需单独启动 `PowerShell` 了。

**note**

> 在 Pycharm 中使用 Powershell 可以点击 Pycharm 底部菜单栏的终端启动 Powershell（小箭头选择）

## 项目开发相关知识

### Git 相关

下面我们将简单介绍 `git` 相关的命令以方便新手学习如何使用 `git` 工具进行项目代码的管理。

对于项目开发来说，经常会涉及到多人合作开发一个项目，因此需要一种工具来管理代码开发的整个过程，涉及代码的修改历史，版本回溯，分支合并等各个方面，这也就是 `git` 工具产生的原因，而 `github` 则是一个代码托管的平台。如果对 `git` 及 `github` 的原理及更多应用感兴趣，可以自行搜索学习。

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

[5. GitHub fork 与 pull request](https://blog.csdn.net/benben0729/article/details/83031135)

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

对于少数没有上传到 conda 的 Python 包，当使用 pip 进行安装时，也需要使用国内的镜像资源进行加速，方法是在用户目录下新建`.pip/pip.conf`文件（先建一个`.pip`文件夹），然后修改如下（使用阿里镜像）：

```
[global]
index-url = https://mirrors.aliyun.com/pypi/simple
```

## 学习资料推荐

### Python 编程语言

#### 1.《Python编程从入门到实践》
<div align="center"><img src="img_16.png"></div>
本书适合完全0基础的读者作为学习Python的入门书籍，全书通过对Python基本语法和简单案例的讲解，能让读者快速掌握Python的核心编程概念，并增强学习Python的信心和兴趣。但是，本书作为0基础的入门书籍，涉及的Python知识点不够，没有着力于培养面向对象编程的思想，读者在理解基础语法后，需要结合[《Python基础教程》](#2Python基础教程)补充知识点与实战训练。此外，书中的实战项目不太实用，可以略过。

>**建议阅读章节：第1章—第10章。**

- [点击下载《Python编程从入门到实践》代码](https://download.ituring.com.cn/book/download/646f556e-220c-4aab-824f-0d2ea7d61a87)

#### 2.《Python基础教程》
<div align="center"><img src="img_17.png"></div>
本书适合对编程语言有过一定接触的读者作为学习Python的入门书籍。全书囊括了Python所有重要的知识点，内容难度适中，并且作者着重强调抽象意识，有助于读者理解面向对象编程，书中也有部分有价值的实战项目，可以用于在了解基础语法之后体会实际的项目编程逻辑和控制。但美中不足的是本书中文版翻译质量较差，建议结合原书对照阅读。

>**建议阅读章节：第1章—第11章，第20章，第22章。**

- [点击下载《Python基础教程》英文版](https://libgen.rocks/get.php?md5=1ad21aa4a661ab666476ee2bc7d499cf&key=9CRH30QJRYL3VQ6Y)

- [点击下载《Python基础教程》代码](https://download.ituring.com.cn/book/download/5ae84b8c-8bb0-48c7-af26-f749e56c761c)

#### 3. 《流畅的Python》
<div align="center"><img src="img_18.png"></div>
Python进阶书籍之一。本书着重讲解Python编程语言独有的特性和设计思路，从Python的角度出发满足各类需求，能够培养读者以Python的风格进行编程。全书分为序章、数据结构、函数、类、控制流以及元编程六个部分，其中前三部分对于数据科学方向的Python实践帮助很大，后三部分更适合用于学习Python软件开发。如果在学习完前2本推荐教材之后直接阅读本书，可能会有一定挑战性，可以先学习一些相对简单的章节，后续再在项目实践的过程中同步学习。
此外，本书的案例相对偏少，建议配合[《Python Cookbook》](#4Python-Cookbook)同步阅读；如果想深入了解数据结构，可以参考[《算法导论》](#3算法导论)的数据结构相关部分。

>**前期建议阅读章节：第1章—第3章，第5章。**

- [点击下载《流畅的Python》代码](https://github.com/fluentpython/example-code/archive/refs/heads/master.zip)

#### 4.《Python Cookbook》
<div align="center"><img src="img_19.png"></div>
Python进阶书籍之一。本书是典型的字典类教材，提供了解决各类实际问题的框架和技巧，但是略微欠缺系统性。推荐读者在基本了解书中涉及哪些内容后，作为`工具书`使用，如果在工作中遇到相似问题，再进行查阅；也可以作为每日读物，每天学习1-2个解决问题的方法。

- [点击下载《Python Cookbook》代码](https://github.com/dabeaz/python-cookbook/archive/refs/heads/master.zip)

### Python 数据分析

#### 1.《Python数据科学手册》
<div align="center"><img src="img_20.png"></div>
本书主要介绍了Python在数据科学领域的几个常用工具: `Jupyter`,`Numpy`,`Pandas`,`Matplotlib`,`Scikit-Learn`。其中`Jupyter`提供日常数据处理的工作环境，`Numpy`用于储存和操作大型数据，`Pandas`处理带标签的大型数据，`Matplotlib`进行数据可视化，`Scikit-Learn`提供传统机器学习的实现。机器学习部分可以略过，[Python 机器学习](#Python-机器学习)部分有更详细教材。

>**建议阅读章节：第1章—第4章。**

- [点击下载《Python数据科学手册》代码](https://github.com/jakevdp/PythonDataScienceHandbook/archive/refs/heads/master.zip)

#### 2.《利用Python进行数据分析》
<div align="center"><img src="img_21.png"></div>
本书大致内容与[《Python数据科学手册》](#1Python数据科学手册)相同，但本书的作者是`Pandas`的主要开发者，如果想从另一个角度了解`Pandas`，可以选择阅读本书。

- [点击下载《利用Python进行数据分析》代码](https://github.com/wesm/pydata-book/archive/refs/heads/3rd-edition.zip)

### Python 网络爬虫

#### 1.《Python3网络爬虫开发实战》
<div align="center"><img src="img_22.png"></div>
本书对网络爬虫从端到端的各项内容都有所涉及，介绍了大量的Python库，对同一个问题也给出了多种解决方法，并且提供了爬虫实战网页，但内容略显冗杂，有堆料的嫌疑，适合作为工具书使用。

>**建议阅读章节：第1章—第7章。**

#### 2.《Python网络爬虫权威指南》
<div align="center"><img src="img_23.png"></div>
本书作为轻量版的爬虫教程，清晰地阐述了爬虫思路和流程，并对可能遇到的问题给出了简要的解决方法，相较于前书，更着重于思路。

- [点击下载《Python网络爬虫权威指南》代码](https://github.com/REMitchell/python-scraping/archive/refs/heads/master.zip)

### Python 机器学习

#### 1.《吴恩达机器学习2022版》(视频)
<div align=center><img width=700 src="img_24.png"/></div>
吴恩达机器学习公开课十分适合初学者入门机器学习，该课程内容简单，涉及面广，并且只需要读者具备基础的数学知识。通过快速学习该课程，读者能够在短期内对各类机器学习基本概念有所了解。需要强调的是，虽然该课程对于模型概念和框架的理解很有帮助，但是实战内容（代码量）很少，需要结合实战材料一起学习。

- [Coursera官网课程链接](https://www.coursera.org/specializations/machine-learning-introduction)

- [B站Up主转载视频链接](https://www.bilibili.com/video/BV1Pa411X76s)

- [课程代码下载](https://github.com/kaieye/2022-Machine-Learning-Specialization/archive/refs/heads/main.zip)

#### 2.《机器学习实战：基于Scikit-Learn、Keras和TensorFlow》
<div align="center"><img src="img_26.png"></div>
本书专注于使用机器学习库进行实际操作，不过多涉及理论模型和实现方法，通过大量的实战代码，手把手教会读者真正从头开始训练机器学习模型的流程，是很适合初学者入门的实战教材。本书分为传统机器学习，即`Scikit-Learn`和深度学习，即`Keras`和`TensorFlow`两部分，由于项目开发更倾向于使用`PyTorch`框架，所以只建议阅读本书的传统机器学习部分。

>**建议阅读章节：第1章—第9章。**

- [点击下载《机器学习实战：基于Scikit-Learn、Keras和TensorFlow》代码](https://github.com/ageron/handson-ml2/archive/refs/heads/master.zip)

#### 3.《深度学习入门：基于Python的理论与实现》
<div align="center"><img src="img_25.png"></div>
本书短小精悍，仅用`Numpy`和Python标准库就完成了简单深度学习框架的从零实现，对深度学习的核心思想与数学模型也有精妙的阐述，让读者对深度学习框架的底层实现有充分了解，十分适合用于深度学习入门。

>**建议阅读章节：全部。**

- [点击下载《深度学习入门：基于Python的理论与实现》代码](https://github.com/oreilly-japan/deep-learning-from-scratch/archive/refs/heads/master.zip)

#### 4.《动手学深度学习》
<div align="center"><img src="img_27.png"></div>
本书包含了简要的深度学习所需数学知识，各类深度学习模型框架，模型优化方法以及开展深度学习所需硬件的介绍。新版教材使用`PyTorch`框架，由浅入深对近年来诞生的不同深度学习框架以及对应效率进行了讲解。在模型介绍之后，针对图片识别和自然语言处理两大应用也附有专项的练习。同时，教材作者还在B站实时直播授课，组织课程竞赛，讲解论文，十分推荐读者结合B站视频同步学习。

>**建议阅读章节：全部。**

- [课程官网](https://zh-v2.d2l.ai/index.html)

- [PyTorch版教材下载](https://zh-v2.d2l.ai/d2l-zh-pytorch.pdf)

- [作者B站主页](https://space.bilibili.com/1567748478)

### 拓展书目

#### 1.《C/C++ Primer Plus》
<div align="center"><img src="img_28.png"></div>
由于Python语言运行效率和内存的限制，以及大量库文件都用C/C++语言实现，在实际工作中，不论是想提升工作代码的效率还是想阅读库文件源码，仅依靠Python是无法满足需求的。相较于Python，C语言更加底层，通过学习C语言，程序员更理解相对底层的机制，并提升程序性能。因此，我们建议读者在工作进入正轨之后再学习C语言，以满足中大型项目或者核心库代码的开发需求。

#### 2.《概率导论》
<div align="center"><img src="img_29.png"></div>
由于机器学习本质上是基于概率模型的，所以在阅读文献或者教材时，不可避免地会遇到大量的统计学概念。因此，简要学习概率论，对于理解各类机器学习模型的工作原理以及阅读文献有重要帮助。

#### 3.《算法导论》
<div align="center"><img src="img_30.png"></div>
编程语言是解决问题的工具，而算法则是解决问题的方法。本书包含了基本的数据结构以及各种算法问题选编，通过对本书的学习，可以很大程度上拓宽读者利用程序解决问题的思路。虽然通读全书的难度曲线十分陡峭，但是我们还是建议读者至少掌握数据结构部分，以帮助读者在面对不同的需求时调用或编写适当的数据结构以简化问题。

#### 4.《统计学习方法》
<div align="center"><img src="img_31.png"></div>
本书主要介绍了各类传统机器学习方法的数学原理，虽然我们十分不建议初学者把时间花在机器学习的数学解释上，但是如果读者有相关基础，并希望能够从更深的理论角度解释选用当前机器学习模型的原因、优劣，或者改进机器学习模型，那么本书是很值得一看的。

#### 5.《深度学习》
<div align="center"><img src="img_32.png"></div>
与[《统计学习方法》](#4统计学习方法)相似，本书主要介绍了深度学习的数学原理，因此，出于同样的理由，我们把本书作为在对深度学习具备一定程度的理解之后的理论层面进阶教材进行推荐。

## 常用库的在线手册
- [Numpy中文文档](https://www.numpy.org.cn/)
- [Pandas](https://pandas.pydata.org/docs/user_guide/index.html#user-guide)
- [Matplotlib](https://matplotlib.org/stable/tutorials/index.html)
- [Scikit-Learn中文文档](https://www.sklearncn.cn/)
- [PyTorch](https://pytorch.org/docs/stable/index.html)
- [RDKit中文文档](http://rdkit.chenzhaoqiang.com/overview.html)

## 贡献者

[Hui Zhou](https://github.com/Rasic2), [DongZhi Li](https://github.com/kealdoom)

## 版权

本教程的版权归 [GongGroup](https://github.com/GongGroup) 所有，更新维护工作暂由 [Hui Zhou](https://github.com/Rasic2) 及 [DongZhi Li](https://github.com/kealdoom) 负责。
