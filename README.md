# 气象可视化平台 (WeatherViz)

这是一款基于 Python 和 Tkinter 的桌面应用程序，用于可视化 NetCDF (.nc) 格式的气象数据，特别针对 ERA5/ERA5-Land 等数据集进行了优化。用户可以选择中国预定义的行政区域（省、市）或自定义坐标边界进行分析和绘图。

**重要提示**: 作者是一名 Python 初学者，该软件主要在 Gemini AI 的协助下编写完成。代码中可能存在不少错漏之处或不够完善的地方，敬请使用者多多包涵！

---

## 主要功能

* **加载 NetCDF 数据**: 自动处理常见的坐标名称 (如 `valid_time` -> `time`)、经度范围 (0-360 -> -180-180) 和纬度顺序问题。
* **区域选择**: 支持选择全国省市级行政区划，或通过经纬度自定义分析区域。
* **多种图表类型**:
    * 时间段平均图 (地图)
    * 风场图 (地图)
    * 时间序列趋势图 (折线图，可选日/月平均聚合)
    * 时间序列箱形图 (按月分布)
    * 日历热力图 (逐日数据可视化)
    * 统计直方图
    * 散点图 (比较两个变量)
* **自定义选项**: 支持自定义图表标题、时间范围选择、地图上显示数值等。
* **输出**: 可预览图表，并选择自动或手动保存为 PNG 图像文件。
* **本地化**: 优先显示中文地名（含港澳台），并支持简繁体转换（需安装 `opencc-python-reimplemented`）。

## 截图

<img width="1744" height="1469" alt="Image" src="https://github.com/user-attachments/assets/4244b25c-2897-441a-af68-58a893c9c788" />

## 运行环境要求

* Python 3.8 或更高版本
* 核心依赖库:
    * `tkinter` (通常 Python 自带)
    * `xarray` (及后端如 `netcdf4` 或 `h5netcdf`)
    * `matplotlib`
    * `cartopy`
    * `geopandas`
    * `regionmask`
    * `numpy`
    * `pandas`
    * `opencc-python-reimplemented` (用于简繁转换和港澳台名称映射)
* 打包依赖 (如果你想自己打包):
    * `pyinstaller`
    * `pyinstaller-hooks-contrib`

## 安装

推荐使用 Conda 创建独立的虚拟环境：

1.  **创建并激活环境** (以 Python 3.9 为例):
    ```bash
    conda create -n weather_map python=3.9
    conda activate weather_map
    ```
2.  **安装主要依赖** (使用 conda-forge 频道):
    ```bash
    conda install -c conda-forge xarray netcdf4 matplotlib cartopy geopandas regionmask numpy pandas pyinstaller pyinstaller-hooks-contrib
    ```
    *(注意: 安装 `cartopy` 和 `geopandas` 可能需要一些时间，因为它们依赖较多底层库)*
3.  **安装 OpenCC**:
    ```bash
    pip install opencc-python-reimplemented
    ```

## 如何使用

### 使用打包好的版本

1.  前往本仓库的 **[Releases 页面](https://github.com/Kaltsit78/weather-visualization-platform)** 下载最新打包好的版本（通常是一个 `.zip` 文件，例如 `WeatherViz_vX.X.X_windows.zip`）。
2.  解压缩下载的 `.zip` 文件到一个你喜欢的目录。
3.  进入解压后的文件夹 (例如 `气象可视化平台_by_Kaltsit`)。
4.  双击运行其中的可执行文件 (例如 `气象可视化平台_by_Kaltsit.exe`)。

## 数据要求

* **气象数据**: 需要 NetCDF (`.nc`) 格式的文件。推荐使用来自[哥白尼气候变化服务中心 (C3S)](https://cds.climate.copernicus.eu/datasets)。数据应至少包含 `time`, `latitude`, `longitude` 这三个维度（或程序可识别的变体如 `valid_time`）。
* **地图数据**: 程序运行所需的中国省界 (`*_1.shp`) 和市界 (`*_2.shp`) Shapefile 文件已包含在项目的 `data` 文件夹中。请确保在运行脚本或打包程序时，`data` 文件夹与主程序文件位于同一目录下。

## 反馈与联系

由于作者是 Python 初学者且主要借助 AI 完成此项目，程序中难免存在错误、考虑不周或性能欠佳之处。

如果您在使用中遇到任何问题、发现 Bug，或有任何改进建议，非常欢迎通过以下方式联系作者 Kaltsit:

* **GitHub Issues**: 在本仓库的 [Issues 页面](https://github.com/Kaltsit78/weather-visualization-platform/issues) 创建一个新的 Issue。 (推荐)
* **Email**: `kaltsit119@gmail.com`

您的反馈对改进这个项目非常有帮助，感谢您的理解与支持！

## 许可证

本项目采用 MIT 许可证授权。详情请见 `LICENSE` 文件。

---
---

# Weather Visualization Platform (WeatherViz)

This is a Python and Tkinter based desktop application for visualizing meteorological data from NetCDF (.nc) files, especially optimized for datasets like ERA5/ERA5-Land. Users can select predefined administrative regions across China (provinces, cities) or define custom coordinate boundaries for analysis and plotting.

**Important Note**: The author is a Python beginner, and this software was primarily developed with the assistance of the Gemini AI. There may be errors, omissions, or imperfections in the code. Your understanding is greatly appreciated!

---

## Features

* **Load NetCDF Data**: Automatically handles common coordinate names (e.g., `valid_time` -> `time`), longitude ranges (0-360 -> -180-180), and latitude order issues.
* **Region Selection**: Supports selecting provincial and city-level administrative regions in China, or defining a custom analysis area using latitude/longitude coordinates.
* **Multiple Chart Types**:
    * Time-averaged maps
    * Wind field maps
    * Time series line charts (with optional daily/monthly aggregation)
    * Time series box plots (monthly distribution)
    * Calendar heatmaps (for daily data visualization)
    * Statistical histograms
    * Scatter plots (for comparing two variables)
* **Customization**: Supports custom chart titles, flexible time period selection, optional display of values on maps, etc.
* **Output**: Allows previewing charts and choosing between automatic or manual saving as PNG image files.
* **Localization**: Prioritizes displaying Chinese place names (including Hong Kong, Macao, Taiwan) and supports Traditional-to-Simplified Chinese conversion (requires `opencc-python-reimplemented`).

## Screenshot

<img width="1744" height="1469" alt="Image" src="https://github.com/user-attachments/assets/4244b25c-2897-441a-af68-58a893c9c788" />

## Requirements

* Python 3.8 or higher
* Core Dependencies:
    * `tkinter` (usually included with Python)
    * `xarray` (with backend like `netcdf4` or `h5netcdf`)
    * `matplotlib`
    * `cartopy`
    * `geopandas`
    * `regionmask`
    * `numpy`
    * `pandas`
    * `opencc-python-reimplemented` (for T<->S conversion and HK/Macao name mapping)
* Packaging Dependencies (if building yourself):
    * `pyinstaller`
    * `pyinstaller-hooks-contrib`

## Installation

Using Conda to create an isolated environment is recommended:

1.  **Create and activate environment** (e.g., with Python 3.9):
    ```bash
    conda create -n weather_map python=3.9
    conda activate weather_map
    ```
2.  **Install main dependencies** (using the conda-forge channel):
    ```bash
    conda install -c conda-forge xarray netcdf4 matplotlib cartopy geopandas regionmask numpy pandas pyinstaller pyinstaller-hooks-contrib
    ```
    *(Note: Installing `cartopy` and `geopandas` might take some time due to their complex dependencies)*
3.  **Install OpenCC**:
    ```bash
    pip install opencc-python-reimplemented
    ```

## Usage

### Use the packaged version

1.  Go to the **[Releases Page](https://github.com/Kaltsit78/weather-visualization-platform)** of this repository and download the latest packaged version (usually a `.zip` file, e.g., `WeatherViz_vX.X.X_windows.zip`).
2.  Extract the downloaded `.zip` file to a location of your choice.
3.  Navigate into the extracted folder (e.g., `气象可视化平台_by_Kaltsit`).
4.  Double-click the executable file (e.g., `气象可视化平台_by_Kaltsit.exe`) to run the application.

## Data Requirements

* **Meteorological Data**: Requires NetCDF (`.nc`) files. Datasets from the [Copernicus Climate Change Service (C3S)](https://cds.climate.copernicus.eu/datasets), such as ERA5 or ERA5-Land, are recommended. The data should include at least the dimensions `time`, `latitude`, and `longitude` (or recognizable variants like `valid_time`).
* **Map Data**: The necessary Shapefile files for China's provincial (`*_1.shp`) and city (`*_2.shp`) boundaries are included in the project's `data` folder. Ensure this `data` folder is present in the same directory as the main script or the packaged executable when running the application.

## Feedback and Contact

As the author is a Python beginner and this project was largely developed with AI assistance, errors, oversights, or suboptimal performance may exist.

If you encounter any issues, find bugs, or have suggestions for improvement, your feedback is highly welcome! Please contact the author, Kaltsit, via:

* **GitHub Issues**: Create a new issue on this repository's [Issues Page](https://github.com/Kaltsit78/weather-visualization-platform/issues). (Preferred)
* **Email**: `kaltsit119@gmail.com`

Your feedback is invaluable for improving this project. Thank you for your understanding and support!

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
