import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import threading
import os
import sys
import glob
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib as mpl
import numpy as np
import pandas as pd

# 尝试导入 geopandas 和 opencc
try:
    import geopandas as gpd
except ImportError:
    gpd = None  # 稍后检查
try:
    from opencc import OpenCC  # <--- 导入 OpenCC
except ImportError:
    OpenCC = None  # 稍后检查

import regionmask
from datetime import datetime
import calendar
import re
import traceback
import subprocess
import configparser
from pathlib import Path

# --- 全局设定：简体中文字体 ---
font_names = ['SimHei', 'Microsoft YaHei', 'sans-serif']
mpl.rcParams['font.sans-serif'] = font_names
mpl.rcParams['axes.unicode_minus'] = False

# --- 创建 OpenCC 转换器实例 ---
if OpenCC:
    try:
        cc = OpenCC('t2s')  # <--- 【修正】: 使用 't2s' 而不是 't2s.json'
    except Exception as e:
        print(f"Error initializing OpenCC: {e}")  # 添加错误处理
        cc = None  # 初始化失败则禁用转换
else:
    cc = None
# --- 结束 ---


# --- 中英对照翻译字典 (用于UI显示) ---
VARIABLE_TRANSLATIONS = {
    "t2m": "2米气温", "tp": "总降水量", "precip": "总降水量", "u10": "10米U风分量",
    "v10": "10米V风分量", "d2m": "2米露点温度", "sp": "地表气压", "msl": "平均海平面气压",
    "e": "蒸发量", "pev": "潜在蒸发量", "ro": "总径流量", "sro": "地表径流量",
    "ssro": "地表下径流量", "sf": "总降雪量",
}

# --- 图表类型中文翻译字典 (用于标题) ---
CHART_TYPE_TRANSLATIONS = {
    "时间段平均图 (地图)": "时间段平均图",
    "风场图": "风场图",
    "时间序列趋势图 (折线)": "时间序列趋势图",
    "时间序列箱形图": "时间序列箱形图",
    "日历热力图": "日历热力图",
    "统计直方图": "统计直方图",
    "散布图": "散布图",
}

# --- 单位格式化字典 ---
UNIT_FORMATTER = {
    "m s**-1": "m/s",
    "kg m**-2 s**-1": "kg/m²/s",
    "m of water equivalent": "m",
    "K": "K",
    "Pa": "Pa",
    "hPa": "hPa",
}

# --- 常见的中文/拼音属性名列表 ---
PROVINCE_CN_ATTRS = ['NL_NAME_1', 'NAME_CHN_1', 'NAME_CN_1', '省名', '省']
PROVINCE_PY_ATTRS = ['NAME_1']
CITY_CN_ATTRS = ['NL_NAME_2', 'NAME_CHN_2', 'NAME_CN_2', '市名', '市', 'NAME_LOC']
CITY_PY_ATTRS = ['NAME_2']
# --- 结束 ---

# --- 港澳台特殊名称映射字典 ---
SPECIAL_AREA_NAME_MAP = {
    "hong kong": "香港特别行政区",
    "macao": "澳门特别行政区",
    "macau": "澳门特别行政区",
    "taiwan": "台湾省",
    # --- 香港下属区域示例 ---
    "central and western": "中西区", "eastern": "东区", "southern": "南区",
    "wan chai": "湾仔区", "kowloon city": "九龙城区", "kwun tong": "观塘区",
    "sham shui po": "深水埗区", "wong tai sin": "黄大仙区", "yau tsim mong": "油尖旺区",
    "islands": "离岛区", "kwai tsing": "葵青区", "north": "北区",
    "sai kung": "西贡区", "sha tin": "沙田区", "tai po": "大埔区",
    "tsuen wan": "荃湾区", "tuen mun": "屯门区", "yuen long": "元朗区",
    # --- 澳门下属区域示例 ---
    "concelho das ilhas": "离岛", "concelho de macau": "澳门半岛",
    "nossa senhora de fatima": "花地玛堂区", "santo antonio": "花王堂区",
    "sao lazaro": "望德堂区", "se": "大堂区", "sao lourenco": "风顺堂区",
    "nossa senhora do carmo": "嘉模堂区", "sao francisco xavier": "圣方济各堂区",
}


# --- 结束 ---


# --- 单位格式化辅助函数 ---
def format_unit_string(unit_str):
    return UNIT_FORMATTER.get(unit_str, unit_str)


# --- 在图表右下角添加数据来源的辅助函数 ---
def add_source_annotation(fig, ds, var_list, dataset_name):
    sources = [ds[var].attrs.get('long_name', var) for var in var_list if var in ds]
    unique_sources = list(dict.fromkeys(sources))
    source_text = f"Data Source: {dataset_name} - {' & '.join(unique_sources)}"
    fig.text(0.98, 0.015, source_text, ha="right", va="bottom", fontsize=8, color='gray')


# --- 通用名称获取函数 - 增加对 "N/A" 的过滤 ---
def get_name(attributes, cn_attrs, py_attrs, unknown_default):
    invalid_placeholders = {"N/A", "NA", ""}
    for attr in cn_attrs:
        name = attributes.get(attr)
        if name and isinstance(name, str):
            name_stripped = name.strip()
            if name_stripped and name_stripped.upper() not in invalid_placeholders:
                return name_stripped
    for attr in py_attrs:
        name = attributes.get(attr)
        if name and isinstance(name, str):
            name_stripped = name.strip()
            if name_stripped and name_stripped.upper() not in invalid_placeholders:
                return name_stripped
    return unknown_default


# --- 结束 ---


# --- 核心绘图函数模块 ---
# ... (所有 plot_... 函数保持不变) ...
def plot_geographical_map(ds_processed, var, title_text, area_bounds, city_shape, province_shape, show_values,
                          dataset_name):
    data = ds_processed[var].mean(dim='time')
    unit = format_unit_string(ds_processed[var].attrs.get('units', ''))
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(area_bounds, crs=ccrs.PlateCarree())
    mesh = data.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap='viridis', add_colorbar=True,
                                cbar_kwargs={'label': unit})
    data.plot.contour(ax=ax, transform=ccrs.PlateCarree(), colors='white', linewidths=0.5)
    if show_values and data.size < 500:
        for i in range(len(data.latitude)):
            for j in range(len(data.longitude)):
                val = data.values[i, j]
                if not np.isnan(val):
                    bgcolor = mesh.cmap(mesh.norm(val))
                    textcolor = 'white' if sum(bgcolor[:3]) < 1.5 else 'black'
                    ax.text(data.longitude[j], data.latitude[i], f'{val:.2f}', ha='center', va='center',
                            color=textcolor, fontsize=8, transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=0.8)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':', linewidth=0.6)
    if province_shape: ax.add_geometries([province_shape], ccrs.PlateCarree(), edgecolor='black', facecolor='none',
                                         linewidth=1.2)
    if city_shape: ax.add_geometries([city_shape], ccrs.PlateCarree(), edgecolor='red', facecolor='none', linewidth=1.8,
                                     linestyle='--')
    gl = ax.gridlines(draw_labels=True, lw=1, color='gray', alpha=0.5, ls='--')
    gl.top_labels = False
    gl.right_labels = False
    ax.set_title(title_text, fontsize=14, pad=20)
    add_source_annotation(fig, ds_processed, [var], dataset_name)
    fig.tight_layout(rect=[0.05, 0.05, 0.95, 0.92])
    return fig


def plot_wind_map(ds_processed, u_var, v_var, title_text, area_bounds, city_shape, province_shape, show_values,
                  dataset_name):
    u_data = ds_processed[u_var].mean(dim='time')
    v_data = ds_processed[v_var].mean(dim='time')
    magnitude = np.sqrt(u_data ** 2 + v_data ** 2)
    skip = max(1, len(u_data.longitude) // 15)
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(area_bounds, crs=ccrs.PlateCarree())
    mesh = magnitude.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap='coolwarm', add_colorbar=True,
                                     cbar_kwargs={'label': 'Wind Speed (m/s)'})
    ax.barbs(u_data.longitude[::skip], u_data.latitude[::skip], u_data.values[::skip, ::skip],
             v_data.values[::skip, ::skip], transform=ccrs.PlateCarree(), length=6)
    if show_values and magnitude.size < 500:
        for i in range(len(magnitude.latitude)):
            for j in range(len(magnitude.longitude)):
                val = magnitude.values[i, j]
                if not np.isnan(val):
                    bgcolor = mesh.cmap(mesh.norm(val))
                    textcolor = 'white' if sum(bgcolor[:3]) < 1.5 else 'black'
                    ax.text(magnitude.longitude[j], magnitude.latitude[i], f'{val:.1f}', ha='center', va='center',
                            color=textcolor, fontsize=8, transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=0.8)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':', linewidth=0.6)
    if province_shape: ax.add_geometries([province_shape], ccrs.PlateCarree(), edgecolor='black', facecolor='none',
                                         linewidth=1.2)
    if city_shape: ax.add_geometries([city_shape], ccrs.PlateCarree(), edgecolor='red', facecolor='none', linewidth=1.8,
                                     linestyle='--')
    gl = ax.gridlines(draw_labels=True, lw=1, color='gray', alpha=0.5, ls='--')
    gl.top_labels = False
    gl.right_labels = False
    ax.set_title(title_text, fontsize=14, pad=20)
    add_source_annotation(fig, ds_processed, [u_var, v_var], dataset_name)
    fig.tight_layout(rect=[0.05, 0.05, 0.95, 0.92])
    return fig


def plot_line_chart(ds_processed, var, title_text, aggregation, dataset_name):
    time_series = ds_processed[var]
    if "latitude" in time_series.dims: time_series = time_series.mean(dim=["latitude", "longitude"])
    if aggregation in {"按日平均": "1D", "按月平均": "1MS"}:
        time_series = time_series.resample(time={"按日平均": "1D", "按月平均": "1MS"}[aggregation]).mean()
    fig, ax = plt.subplots(figsize=(12, 6))
    time_series.to_series().plot(marker='o', ax=ax)
    unit = format_unit_string(ds_processed[var].attrs.get('units', ''))
    ax.set_title(title_text, fontsize=14)
    ax.set_ylabel(f"Value ({unit})")
    ax.set_xlabel("Time")
    ax.grid(True)
    add_source_annotation(fig, ds_processed, [var], dataset_name)
    fig.tight_layout(rect=[0.05, 0.08, 0.95, 0.92])
    return fig


def plot_boxplot(ds_processed, var, title_text, dataset_name):
    time_series = ds_processed[var]
    if "latitude" in time_series.dims: time_series = time_series.mean(dim=["latitude", "longitude"])
    df = time_series.to_series().to_frame(name=var)
    df['month'] = df.index.month
    fig, ax = plt.subplots(figsize=(12, 6))
    df.boxplot(column=var, by='month', ax=ax)
    unit = format_unit_string(ds_processed[var].attrs.get('units', ''))
    ax.set_title(title_text, fontsize=14)
    fig.suptitle('')
    ax.set_ylabel(f"Daily Value ({unit})")
    ax.set_xlabel("Month")
    ax.grid(True)
    add_source_annotation(fig, ds_processed, [var], dataset_name)
    fig.tight_layout(rect=[0.05, 0.08, 0.95, 0.92])
    return fig


def plot_calendar_heatmap(ds_processed, var, title_text, dataset_name):
    year = ds_processed.time.dt.year.values[0]
    if ds_processed.time.size == 0: raise ValueError(f"找不到 {year}年 的数据！")

    time_series = ds_processed[var]
    if "latitude" in time_series.dims: time_series = time_series.mean(dim=["latitude", "longitude"])

    aggregation_method = "sum" if var in ["tp", "precip", "ro", "sro", "ssro", "e", "pev", "sf"] else "mean"

    try:
        if aggregation_method == "sum":
            year_series = time_series.resample(time="1D", label="left", closed="left").sum().to_series()
        else:
            year_series = time_series.resample(time="1D", label="left", closed="left").mean().to_series()
    except Exception as e:
        raise ValueError(f"数据重采样为'1D'失败: {e}。请检查原始数据的时间坐标。")

    year_series = year_series[year_series.index.year == int(year)]

    if year_series.empty or year_series.isnull().all():
        raise ValueError(f"重采样 {year}年 数据后，所有值都为空。")

    vmin = year_series.min()
    vmax = year_series.max()
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = 'viridis'

    fig = plt.figure(figsize=(15, 13))
    fig.suptitle(title_text, fontsize=14, y=0.98)

    gs_months = fig.add_gridspec(4, 3, left=0.05, right=0.95, top=0.92, bottom=0.1, hspace=0.4, wspace=0.1)

    cbar_ax = fig.add_axes([0.3, 0.05, 0.4, 0.015])
    axes = [fig.add_subplot(gs_months[i, j]) for i in range(4) for j in range(3)]

    for month in range(1, 13):
        ax = axes[month - 1]
        cal = calendar.monthcalendar(int(year), month)
        month_data = np.full((len(cal), 7), np.nan)

        for week_idx, week in enumerate(cal):
            for day_idx, day in enumerate(week):
                if day != 0:
                    try:
                        current_date = pd.to_datetime(f'{year}-{month}-{day}')
                        val = year_series.get(current_date)
                        if pd.notnull(val):
                            month_data[week_idx, day_idx] = val
                    except (KeyError, ValueError, IndexError):
                        pass

        im = ax.imshow(month_data, cmap=cmap, norm=norm, interpolation='none', aspect='equal')

        ax.set_title(f'{month}月', fontsize=10)
        ax.set_xticks(np.arange(7))
        ax.set_xticklabels(['一', '二', '三', '四', '五', '六', '日'], fontsize=8)
        ax.set_yticks(np.arange(len(cal)))
        ax.set_yticklabels([])
        ax.tick_params(axis='y', which='both', left=False, right=False)
        ax.tick_params(axis='x', which='both', bottom=False, top=False)

    unit = format_unit_string(ds_processed[var].attrs.get('units', ''))
    if aggregation_method == "sum" and unit:
        unit = f"{unit}/day"

    mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(mappable, cax=cbar_ax, orientation='horizontal')
    cbar.set_label(f"逐日{aggregation_method.title()} ({unit})", fontsize=8)
    cbar.ax.tick_params(labelsize=6)

    add_source_annotation(fig, ds_processed, [var], dataset_name)

    return fig


def plot_histogram(ds_processed, var, title_text, dataset_name):
    data_flat = ds_processed[var].values.flatten()
    data_flat = data_flat[~np.isnan(data_flat)]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(data_flat, bins=50, density=True)
    unit = format_unit_string(ds_processed[var].attrs.get('units', ''))
    ax.set_title(title_text, fontsize=14)
    ax.set_xlabel(f"Value ({unit})")
    ax.set_ylabel("Density")
    ax.grid(True)
    add_source_annotation(fig, ds_processed, [var], dataset_name)
    fig.tight_layout(rect=[0.05, 0.08, 0.95, 0.9])
    return fig


def plot_scatter(ds_processed, var1, var2, title_text, dataset_name):
    var1_series = ds_processed[var1]
    if "latitude" in var1_series.dims: var1_series = var1_series.mean(dim=["latitude", "longitude"])
    var2_series = ds_processed[var2]
    if "latitude" in var2_series.dims: var2_series = var2_series.mean(dim=["latitude", "longitude"])
    unit1 = format_unit_string(ds_processed[var1].attrs.get('units', ''))
    unit2 = format_unit_string(ds_processed[var2].attrs.get('units', ''))
    long_name1 = ds_processed[var1].attrs.get('long_name', var1)
    long_name2 = ds_processed[var2].attrs.get('long_name', var2)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(var1_series, var2_series, alpha=0.5)
    ax.set_title(title_text, fontsize=14)
    ax.set_xlabel(f"{long_name1} ({unit1})")
    ax.set_ylabel(f"{long_name2} ({unit2})")
    ax.grid(True)
    add_source_annotation(fig, ds_processed, [var1, var2], dataset_name)
    fig.tight_layout(rect=[0.05, 0.08, 0.95, 0.9])
    return fig


# --- 图形使用者介面 (GUI) ---
class LocalWeatherPlatform:
    CONFIG_FILE = Path.home() / '.weather_app_config.ini'

    # --- 替换现有的 __init__ 方法 ---
    def __init__(self, root):
        self.root = root
        self.root.title("气象可视化平台 V17.26")  # 保持或更新版本号

        # --- 设置初始窗口尺寸 ---
        window_width_init = 1400
        window_height_init = 1150  # 使用新高度
        self.root.geometry(f"{window_width_init}x{window_height_init}")
        # --- 结束 ---

        # --- 【新】更可靠的窗口居中 ---
        self.root.withdraw()  # 1. 先隐藏窗口
        self.root.update_idletasks()  # 2. 处理初始尺寸请求

        # 3. 计算居中位置
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        center_x = int(screen_width / 2 - window_width_init / 2)
        # 向上偏移以避开任务栏
        offset_y = 45  # 可以调整这个值 (例如 40, 50, 60)
        center_y = int(screen_height / 2 - window_height_init / 2) - offset_y
        center_y = max(0, center_y)  # 确保不移出屏幕顶部

        # 4. 设置最终尺寸和计算出的位置
        self.root.geometry(f'{window_width_init}x{window_height_init}+{center_x}+{center_y}')

        # 5. 显示窗口
        self.root.deiconify()
        # --- 【新】结束 ---

        self.ds = None
        self.city_shapes = {}
        self.province_shapes = {}
        self.city_names = ["--- 正在加载地图数据 ---"]

        self.last_year = tk.StringVar()
        self.last_month = tk.StringVar()
        self.last_day = tk.StringVar()
        self.last_hour = tk.StringVar()
        self.last_city = tk.StringVar(value="自定义区域")
        self.last_coords = tk.StringVar(value="35, 105, 20, 125")

        self.save_mode_var = tk.StringVar()
        self._load_config()

        self.current_fig = None
        self.current_output_filename = None

        main_frame = ttk.Frame(root, padding="10")
        main_frame.pack(fill="both", expand=True)

        bottom_frame = ttk.Frame(main_frame)
        bottom_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=(5, 0))

        author_label = ttk.Label(bottom_frame, text="作者: Kaltsit (kaltsit119@gmail.com)", font=("Calibri", 9))
        author_label.pack(side=tk.BOTTOM, fill="x")

        self.status_var = tk.StringVar(value="正在背景初始化地图数据，请稍候...")
        status_label = ttk.Label(bottom_frame, textvariable=self.status_var,
                                 wraplength=window_width_init - 40)  # 使用初始宽度计算
        status_label.pack(side=tk.BOTTOM, fill="x")

        paned_window = ttk.PanedWindow(main_frame, orient=tk.HORIZONTAL)
        paned_window.pack(fill="both", expand=True, pady=(0, 5))

        left_controls_frame = ttk.Frame(paned_window, padding="5")

        load_frame = ttk.LabelFrame(left_controls_frame, text="1. 加载天气数据")
        load_frame.pack(fill="x", pady=5)
        self.file_path_var = tk.StringVar()
        ttk.Entry(load_frame, textvariable=self.file_path_var, width=50).pack(side="left", fill="x", expand=True,
                                                                              padx=5, pady=5)
        self.info_button = ttk.Button(load_frame, text="(!)", command=self.show_data_instructions, width=3)
        self.info_button.pack(side="right", padx=(0, 5), pady=5)
        self.load_nc_button = ttk.Button(load_frame, text="加载数据", command=self.load_netcdf, state="disabled")
        self.load_nc_button.pack(side="right", padx=5, pady=5)

        self.settings_frame = ttk.LabelFrame(left_controls_frame, text="2. 参数设定")
        self.settings_frame.pack(fill="x", pady=5)
        ttk.Label(self.settings_frame, text="数据集名称:").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.dataset_name_var = tk.StringVar(value="ERA5-Land")
        ttk.Entry(self.settings_frame, textvariable=self.dataset_name_var).grid(row=0, column=1, columnspan=3,
                                                                                sticky="ew", padx=5, pady=2)
        ttk.Label(self.settings_frame, text="图表类型:").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        self.chart_type_var = tk.StringVar()
        self.chart_type_cb = ttk.Combobox(self.settings_frame, textvariable=self.chart_type_var, state="disabled",
                                          values=list(CHART_TYPE_TRANSLATIONS.keys()))
        self.chart_type_cb.grid(row=1, column=1, columnspan=3, sticky="ew", padx=5, pady=2)
        self.chart_type_cb.bind("<<ComboboxSelected>>", self.toggle_options)

        ttk.Label(self.settings_frame, text="自定义标题:").grid(row=2, column=0, sticky="w", padx=5, pady=2)
        self.custom_title_var = tk.StringVar()
        ttk.Entry(self.settings_frame, textvariable=self.custom_title_var).grid(row=2, column=1, columnspan=2,
                                                                                sticky="ew", padx=5, pady=2)
        ttk.Label(self.settings_frame, text="(留空则自动生成)", foreground="gray").grid(row=2, column=3, sticky="w",
                                                                                        padx=5)

        self.var_frame = ttk.Frame(self.settings_frame)
        self.var_frame.grid(row=3, column=0, columnspan=4, sticky="ew", pady=2)
        self.agg_frame = ttk.Frame(self.settings_frame)
        self.agg_frame.grid(row=4, column=0, columnspan=4, sticky="ew", pady=2)
        ttk.Label(self.settings_frame, text="分析区域:").grid(row=5, column=0, sticky="w", padx=5, pady=2)
        self.city_var = self.last_city
        self.city_cb = ttk.Combobox(self.settings_frame, textvariable=self.city_var, values=self.city_names,
                                    state="disabled")
        self.city_cb.grid(row=5, column=1, columnspan=3, sticky="ew", padx=5, pady=2)
        self.city_cb.bind("<<ComboboxSelected>>", self.toggle_area_input)
        self.area_entry_label = ttk.Label(self.settings_frame, text="自定义坐标:")
        self.area_var = self.last_coords
        self.area_entry = ttk.Entry(self.settings_frame, textvariable=self.area_var, state="disabled")
        self.area_info_label = ttk.Label(self.settings_frame, text="(北, 西, 南, 东)")
        self.time_frame = ttk.Frame(self.settings_frame)
        self.time_frame.grid(row=7, column=0, columnspan=4, sticky="ew", pady=2)
        self.show_values_var = tk.BooleanVar(value=False)
        self.show_values_check = ttk.Checkbutton(self.settings_frame, text="在地图上显示数值",
                                                 variable=self.show_values_var)
        self.show_values_check.grid(row=8, column=0, columnspan=2, padx=5, pady=5, sticky='w')

        self.settings_frame.grid_columnconfigure(1, weight=1)
        self.settings_frame.grid_columnconfigure(2, weight=1)

        action_frame = ttk.LabelFrame(left_controls_frame, text="3. 执行与保存")
        action_frame.pack(fill="x", pady=10)

        path_row = ttk.Frame(action_frame)
        path_row.pack(fill='x', expand=True, pady=(0, 5))

        ttk.Label(path_row, text="保存路径:").grid(row=0, column=0, sticky='w', padx=5, pady=2)
        self.save_path_var = tk.StringVar(value=self.last_save_path)

        path_entry = ttk.Entry(path_row, textvariable=self.save_path_var, state='normal')
        path_entry.grid(row=0, column=1, sticky='ew', padx=5, pady=2)
        path_entry.bind("<FocusOut>", self._on_path_entry_change)
        path_entry.bind("<Return>", self._on_path_entry_change)

        ttk.Button(path_row, text="选择路径", command=self.select_save_path, width=10).grid(row=0, column=2, sticky='e',
                                                                                            padx=2, pady=2)
        ttk.Button(path_row, text="打开文件夹", command=self._open_save_folder, width=10).grid(row=0, column=3,
                                                                                               sticky='e', padx=2,
                                                                                               pady=2)

        path_row.grid_columnconfigure(1, weight=1)
        path_row.grid_columnconfigure(0, weight=0)
        path_row.grid_columnconfigure(2, weight=0)
        path_row.grid_columnconfigure(3, weight=0)

        name_row = ttk.Frame(action_frame)
        name_row.pack(fill='x', expand=True, pady=(0, 5))

        ttk.Label(name_row, text="自定义文件名:").grid(row=0, column=0, sticky='w', padx=5, pady=2)
        self.custom_filename_var = tk.StringVar()
        ttk.Entry(name_row, textvariable=self.custom_filename_var).grid(row=0, column=1, sticky='ew', padx=5, pady=2)

        ttk.Label(name_row, text="(.png, 留空则默认)", foreground="gray").grid(row=0, column=2, sticky='w', padx=5,
                                                                               pady=2)

        name_row.grid_columnconfigure(1, weight=1)
        name_row.grid_columnconfigure(0, weight=0)
        name_row.grid_columnconfigure(2, weight=0)

        mode_frame = ttk.Frame(action_frame)
        mode_frame.pack(fill='x', pady=(5, 0))
        ttk.Radiobutton(mode_frame, text="预览并自动保存", variable=self.save_mode_var, value="auto",
                        command=self._save_config).pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(mode_frame, text="预览后手动保存", variable=self.save_mode_var, value="manual",
                        command=self._save_config).pack(side=tk.LEFT, padx=5)

        button_frame = ttk.Frame(action_frame)
        button_frame.pack(fill='x', pady=(5, 5))

        self.generate_button = ttk.Button(button_frame, text="生成图表", command=self.start_task, state="disabled")
        self.generate_button.grid(row=0, column=0, sticky='ew', ipady=5, padx=(5, 2))

        self.save_chart_button = ttk.Button(button_frame, text="保存图表", command=self._save_current_chart,
                                            state="disabled")
        self.save_chart_button.grid(row=0, column=1, sticky='ew', ipady=5, padx=(2, 5))

        button_frame.grid_columnconfigure(0, weight=2)
        button_frame.grid_columnconfigure(1, weight=1)

        self.preview_frame = ttk.LabelFrame(paned_window, text="4. 图表预览")

        # 保持宽度权重不变 (1:5)
        paned_window.add(left_controls_frame, weight=1)
        paned_window.add(self.preview_frame, weight=5)

        self.placeholder_label = ttk.Label(self.preview_frame, text="图表将在此处预览")
        self.placeholder_label.pack(pady=50)
        self.canvas_widget = None
        self.toolbar_widget = None

        threading.Thread(target=self.load_shapefiles_startup, daemon=True).start()

    # --- __init__ 方法结束 ---

    # ... (所有其他方法保持不变) ...

    def _load_config(self):
        config = configparser.ConfigParser()
        self.last_save_path = os.getcwd()  # Default
        save_mode_default = 'auto'  # Default save mode
        try:
            if self.CONFIG_FILE.exists():
                config.read(self.CONFIG_FILE, encoding='utf-8')
                path = config.get('Settings', 'save_path', fallback=os.getcwd())
                self.last_save_path = path
                save_mode_default = config.get('Settings', 'save_mode', fallback='auto')
        except Exception as e:
            print(f"读取配置文件失败: {e}")
        self.save_mode_var.set(save_mode_default)

    def _save_config(self):
        config = configparser.ConfigParser()
        if not config.has_section('Settings'): config.add_section('Settings')
        config.set('Settings', 'save_path', self.save_path_var.get())
        config.set('Settings', 'save_mode', self.save_mode_var.get())
        try:
            with open(self.CONFIG_FILE, 'w', encoding='utf-8') as configfile:
                config.write(configfile)
        except Exception as e:
            print(f"保存配置文件失败: {e}")

    def _on_path_entry_change(self, event=None):
        self._save_config()
        new_path = self.save_path_var.get()
        if not os.path.isdir(new_path):
            self.status_var.set(f"提示: 您输入的路径 '{new_path}' 当前不存在。")
        else:
            self.status_var.set(f"保存路径已更新为: {new_path}")

    def select_save_path(self):
        path = filedialog.askdirectory(initialdir=self.save_path_var.get())
        if path:
            self.save_path_var.set(path)
            self._save_config()

    def _open_save_folder(self):
        folder_path = os.path.abspath(self.save_path_var.get())
        if not os.path.isdir(folder_path):
            messagebox.showwarning("路径无效", f"找不到文件夹:\n{folder_path}")
            return
        try:
            if sys.platform == "win32":
                os.startfile(folder_path)
            elif sys.platform == "darwin":
                subprocess.Popen(["open", folder_path])
            else:
                subprocess.Popen(["xdg-open", folder_path])
        except Exception as e:
            messagebox.showerror("打开失败", f"无法打开文件夹：\n{e}")

    def show_data_instructions(self):
        title = "数据下载注意事项"
        message = """为了确保本程序能正确分析您的数据，请在下载时注意以下几点：

● 数据源推荐:
  - 推荐从「哥白尼气候变化服务中心(C3S)」网站下载。
  - 核心数据集推荐「ERA5」或「ERA5-Land」。

● 分辨率选择:
  - 若要进行城市级别的精细分析，请选择高分辨率的数据，例如 0.1° 或 0.25°。
  - 若使用低分辨率数据(如1°)，请选择较大的分析区域(如整个省或自定义大范围坐标)。

● 【【【 最重要的一点 】】】
  - 在下载网站上，文件格式(Format)选项务必选择 NetCDF (.nc)！

● 本程序很智能！
  - 程序会自动识别并修复不规范的数据格式，包括：
    - 坐标名称不符 (如 valid_time 会自动转为 time)
    - 经度范围不符 (0-360° 会自动转为 -180-180°)
    - 纬度顺序颠倒 (从北到南会自动转为从南到北)

--------------------------------------------------
作者: Kaltsit (kaltsit119@gmail.com)
如果您在使用中遇到任何问题或有改进建议，欢迎反馈至我的邮箱！
"""
        messagebox.showinfo(title, message)

    def resource_path(self, relative_path):
        try:
            base_path = sys._MEIPASS
        except Exception:
            base_path = os.path.abspath(".")
        return os.path.join(base_path, relative_path)

    # --- 【最终修正】: 显式添加省份 + 映射 + 清理 + 过滤 + 繁转简 ---
    def load_shapefiles_startup(self):
        if not cc:
            self.root.after(0, lambda: messagebox.showwarning("缺少库",
                                                              "未找到 'opencc-python-reimplemented' 库。\n繁体字地区名称将无法自动转换为简体。\n请运行 'pip install opencc-python-reimplemented' 安装。"))

        try:
            shape_dir = self.resource_path('data')
            province_files = glob.glob(os.path.join(shape_dir, '*_1.shp'))
            city_files = glob.glob(os.path.join(shape_dir, '*_2.shp'))
            if not city_files or not province_files: raise ValueError(
                "在软体内部找不到省份 (*_1.shp) 或城市 (*_2.shp) 的地图文件。")

            all_region_data = {}
            province_geometries = {}

            # --- 1. 处理省份文件 (*_1.shp) ---
            for shp in province_files:
                reader = shpreader.Reader(shp, encoding='utf-8')
                for rec in reader.records():
                    province_name_raw = get_name(rec.attributes, PROVINCE_CN_ATTRS, PROVINCE_PY_ATTRS, "未知省份")
                    province_name_cleaned = province_name_raw.split('|')[-1].strip()
                    province_name_mapped = SPECIAL_AREA_NAME_MAP.get(province_name_cleaned.lower(),
                                                                     province_name_cleaned)
                    province_name_simplified = cc.convert(province_name_mapped) if cc else province_name_mapped  # 繁转简

                    if province_name_simplified != "未知省份":
                        display_name = province_name_simplified
                        province_geometries[province_name_cleaned] = rec.geometry
                        all_region_data[display_name] = {
                            'name': display_name,
                            'shape': rec.geometry,
                            'province': province_name_cleaned
                        }
                    else:
                        print(f"Skipping unknown province in _1.shp: {rec.attributes}")

            # --- 2. 处理城市文件 (*_2.shp) ---
            for shp in city_files:
                reader = shpreader.Reader(shp, encoding='utf-8')
                geoms = list(reader.geometries())
                recs = list(reader.records())
                for i, rec in enumerate(recs):
                    province_name_raw = get_name(rec.attributes, PROVINCE_CN_ATTRS, PROVINCE_PY_ATTRS, "未知省份")
                    city_name_raw = get_name(rec.attributes, CITY_CN_ATTRS, CITY_PY_ATTRS, "未知城市")

                    province_name_cleaned = province_name_raw.split('|')[-1].strip()
                    city_name_cleaned = city_name_raw.split('|')[-1].strip()

                    province_name_mapped = SPECIAL_AREA_NAME_MAP.get(province_name_cleaned.lower(),
                                                                     province_name_cleaned)
                    city_name_mapped = SPECIAL_AREA_NAME_MAP.get(city_name_cleaned.lower(), city_name_cleaned)

                    province_name_simplified = cc.convert(province_name_mapped) if cc else province_name_mapped
                    city_name_simplified = cc.convert(city_name_mapped) if cc else city_name_mapped

                    # --- 过滤掉完全未知的条目 ---
                    if province_name_simplified == "未知省份" and city_name_simplified == "未知城市":
                        print(f"Skipping record with unknown simplified province and city in _2.shp: {rec.attributes}")
                        continue
                    # --- 过滤掉只有省名没有市名的条目 ---
                    if city_name_simplified == "未知城市" and province_name_simplified != "未知省份":
                        continue

                    # --- 格式化显示名称 (使用简体) ---
                    if province_name_simplified == city_name_simplified or province_name_simplified == "未知省份":
                        display_name = city_name_simplified
                    else:
                        display_name = f"{province_name_simplified} / {city_name_simplified}"

                    # 添加到总区域数据中
                    all_region_data[display_name] = {
                        'name': display_name,
                        'shape': geoms[i],
                        'province': province_name_cleaned  # 内部关联键仍然使用清理后的原始省名
                    }

            # --- 3. 整理并排序 ---
            sorted_display_names = sorted(all_region_data.keys())
            self.city_names = ["自定义区域"] + sorted_display_names
            self.city_shapes = all_region_data
            self.province_shapes = province_geometries

            self.root.after(0, self.update_ui_after_shapefile_load)
        except Exception as e:
            traceback.print_exc()
            self.root.after(0, lambda exc=e: messagebox.showerror("地图数据初始化错误",
                                                                  f"软体启动时加载地图文件失败：\n{exc}"))
            self.root.after(0, lambda: self.status_var.set("地图数据初始化失败！"))

    # --- 【最终修正结束】 ---

    def update_ui_after_shapefile_load(self):
        self.city_cb.config(values=self.city_names, state="readonly")
        last_city_value = self.last_city.get()
        if last_city_value in self.city_names:
            self.city_cb.set(last_city_value)
        else:
            self.city_cb.set(self.city_names[0])
            self.last_city.set(self.city_names[0])

        self.toggle_area_input()
        self.load_nc_button.config(state="normal")
        self.status_var.set(f"地图数据初始化成功！共加载 {len(self.city_names) - 1} 个地区。请加载天气数据。")

    def load_netcdf(self):
        filepath = filedialog.askopenfilename(filetypes=[("NetCDF files", "*.nc")])
        if not filepath: return
        try:
            self.status_var.set(f"正在加载和标准化数据...")
            self.file_path_var.set(filepath)
            ds = xr.open_dataset(filepath, decode_times=True)
            if 'valid_time' in ds.coords: ds = ds.rename({'valid_time': 'time'})
            if 'longitude' in ds.coords and ds.coords['longitude'].max() > 180:
                ds.coords['longitude'] = (ds.coords['longitude'] + 180) % 360 - 180
                ds = ds.sortby(ds.longitude)
            if 'latitude' in ds.coords and ds.coords['latitude'].values[0] > ds.coords['latitude'].values[-1]:
                ds = ds.sortby('latitude')
            self.ds = ds
            self.data_vars_raw = [var for var in self.ds.data_vars if
                                  'time' in self.ds[var].dims and 'latitude' in self.ds[var].dims and 'longitude' in
                                  self.ds[var].dims]
            if not self.data_vars_raw: raise ValueError(
                "标准化后，文件中仍未找到有效的、包含[time, latitude, longitude]维度的气象变数。")
            self.data_vars_display = [f"{VARIABLE_TRANSLATIONS.get(var, var)} ({var})" for var in self.data_vars_raw]
            self.chart_type_cb.config(state="readonly")
            self.chart_type_cb.set(self.chart_type_cb['values'][0])
            self.toggle_options()
            self.generate_button.config(state="normal")
            self.status_var.set(f"数据加载和标准化成功！")
        except Exception as e:
            messagebox.showerror("文件读取错误", f"无法读取或处理 NetCDF 文件：\n{e}")

    def toggle_area_input(self, event=None):
        self.last_city.set(self.city_var.get())
        if self.city_var.get() == "自定义区域": self.last_coords.set(self.area_var.get())
        is_custom = self.city_var.get() == "自定义区域"
        self.area_entry_label.grid(row=6, column=0, sticky="w", padx=5, pady=2)
        self.area_entry.grid(row=6, column=1, columnspan=2, sticky="ew", padx=5, pady=2)
        self.area_info_label.grid(row=6, column=3, sticky="w", padx=5, pady=2)
        self.area_entry.config(state="normal" if is_custom else "disabled")

    def toggle_options(self, event=None):
        if hasattr(self, 'year_cb'): self.last_year.set(self.year_var.get())
        if hasattr(self, 'month_cb'): self.last_month.set(self.month_var.get())
        if hasattr(self, 'day_cb'): self.last_day.set(self.day_var.get())
        if hasattr(self, 'hour_cb'): self.last_hour.set(self.hour_var.get())

        for widget in self.var_frame.winfo_children(): widget.destroy()
        for widget in self.agg_frame.winfo_children(): widget.destroy()
        for widget in self.time_frame.winfo_children(): widget.destroy()

        chart_type = self.chart_type_var.get()
        display_list = self.data_vars_display if hasattr(self, 'data_vars_display') else []
        if chart_type in ["时间段平均图 (地图)", "时间序列趋势图 (折线)", "时间序列箱形图", "日历热力图", "统计直方图"]:
            ttk.Label(self.var_frame, text="气象要素:").pack(side="left", padx=5)
            self.var1_cb = ttk.Combobox(self.var_frame, values=display_list, state="readonly")
            self.var1_cb.pack(side="left", fill="x", expand=True, padx=5)
            if display_list: self.var1_cb.set(display_list[0])
        elif chart_type in ["风场图", "散布图"]:
            text1, text2 = ("U风分量 (东西风):", "V风分量 (南北风):") if chart_type == "风场图" else (
                "变数 X (气象要素):", "变数 Y (气象要素):")
            ttk.Label(self.var_frame, text=text1).pack(side="left", padx=5)
            self.var1_cb = ttk.Combobox(self.var_frame, values=display_list, state="readonly", width=20)
            self.var1_cb.pack(side="left", padx=5)
            ttk.Label(self.var_frame, text=text2).pack(side="left", padx=5)
            self.var2_cb = ttk.Combobox(self.var_frame, values=display_list, state="readonly", width=20)
            self.var2_cb.pack(side="left", padx=5)
            if display_list and len(display_list) > 1: self.var1_cb.set(display_list[0]); self.var2_cb.set(
                display_list[1])
        if chart_type == "时间序列趋势图 (折线)":
            ttk.Label(self.agg_frame, text="时间聚合方式:").pack(side="left", padx=5)
            self.agg_var = tk.StringVar(value="原始频率")
            self.agg_cb = ttk.Combobox(self.agg_frame, textvariable=self.agg_var,
                                       values=["原始频率", "按日平均", "按月平均"], state="readonly")
            self.agg_cb.pack(side="left", padx=5)

        self.year_var = tk.StringVar(value=self.last_year.get())
        self.month_var = tk.StringVar(value=self.last_month.get())
        self.day_var = tk.StringVar(value=self.last_day.get())
        self.hour_var = tk.StringVar(value=self.last_hour.get())

        ttk.Label(self.time_frame, text="时间(选填):").pack(side="left", padx=5)
        self.year_cb = ttk.Combobox(self.time_frame, textvariable=self.year_var,
                                    values=[""] + [str(y) for y in range(1940, datetime.now().year + 3)], width=6)
        self.year_cb.pack(side="left", padx=2)
        ttk.Label(self.time_frame, text="年").pack(side="left")
        if chart_type in ["时间段平均图 (地图)", "风场图", "时间序列趋势图 (折线)"]:
            self.month_cb = ttk.Combobox(self.time_frame, textvariable=self.month_var,
                                         values=[""] + [str(m) for m in range(1, 13)], width=4)
            self.month_cb.pack(side="left", padx=2)
            ttk.Label(self.time_frame, text="月").pack(side="left")
            self.day_cb = ttk.Combobox(self.time_frame, textvariable=self.day_var,
                                       values=[""] + [str(d) for d in range(1, 32)], width=4)
            self.day_cb.pack(side="left", padx=2)
            ttk.Label(self.time_frame, text="日").pack(side="left")
        if chart_type in ["时间段平均图 (地图)", "风场图"]:
            self.hour_cb = ttk.Combobox(self.time_frame, textvariable=self.hour_var,
                                        values=[""] + [str(h) for h in range(24)], width=4)
            self.hour_cb.pack(side="left", padx=2)
            ttk.Label(self.time_frame, text="时").pack(side="left")
        if chart_type in ["时间段平均图 (地图)", "风场图"]:
            self.show_values_check.config(state="normal")
        else:
            self.show_values_check.config(state="disabled");
            self.show_values_var.set(False)

    def get_var_from_display(self, display_string):
        match = re.search(r'\((.*?)\)', display_string)
        return match.group(1) if match else display_string

    def start_task(self):
        if self.current_fig:
            plt.close(self.current_fig)
            self.current_fig = None
            self.current_output_filename = None

        if hasattr(self, 'year_cb'): self.last_year.set(self.year_var.get())
        if hasattr(self, 'month_cb'): self.last_month.set(self.month_var.get())
        if hasattr(self, 'day_cb'): self.last_day.set(self.day_var.get())
        if hasattr(self, 'hour_cb'): self.last_hour.set(self.hour_var.get())
        self.last_city.set(self.city_var.get())
        if self.city_var.get() == "自定义区域": self.last_coords.set(self.area_var.get())

        self.generate_button.config(state="disabled")
        self.save_chart_button.config(state="disabled")

        try:
            params = {
                'region_name': self.city_var.get(),
                'chart_type': self.chart_type_var.get(),
                'dataset_name': self.dataset_name_var.get(),
                'save_path': self.save_path_var.get(),
                'custom_filename': self.custom_filename_var.get(),
                'custom_title': self.custom_title_var.get(),
                'save_mode': self.save_mode_var.get(),
                'area_coords_str': self.area_var.get(),
                'year': self.year_var.get(),
                'month': self.month_var.get(),
                'day': self.day_var.get(),
                'hour': self.hour_var.get(),
                'var1_display': self.var1_cb.get(),
                'var1': self.get_var_from_display(self.var1_cb.get()),
                'show_values': self.show_values_var.get()
            }

            if params['chart_type'] in ["风场图", "散布图"]:
                params['var2_display'] = self.var2_cb.get()
                params['var2'] = self.get_var_from_display(self.var2_cb.get())
            if params['chart_type'] == "时间序列趋势图 (折线)":
                params['aggregation'] = self.agg_var.get()

        except Exception as e:
            messagebox.showerror("参数错误", f"获取界面参数时出错: {e}")
            self.generate_button.config(state="normal")
            self.save_chart_button.config(state="disabled")
            return

        threading.Thread(target=self.run_generation, args=(params,), daemon=True).start()

    def _update_status_safe(self, message):
        self.status_var.set(message)

    def run_generation(self, params):
        try:
            self.root.after(0, self._update_status_safe, "任务开始，正在后台准备数据...")

            if params['region_name'] == "自定义区域":
                coords_str = params['area_coords_str'].split(',')
                north, west, south, east = [float(c.strip()) for c in coords_str]
                if not (north > south and east > west):
                    raise ValueError(f"自定义坐标无效: 北({north})必须大于南({south})，东({east})必须大于西({west})")
                ds_processed = self.ds.sel(latitude=slice(south, north), longitude=slice(west, east))
                params['base_title_part'] = f"自定义区域 ({west:.1f}-{east:.1f}, {south:.1f}-{north:.1f})"
                lon_buffer = (east - west) * 0.1
                lat_buffer = (north - south) * 0.1
                params['area_bounds'] = [west - lon_buffer, east + lon_buffer, south - lat_buffer, north + lat_buffer]
                params['city_shape'], params['province_shape'] = None, None
            else:
                if params['region_name'] not in self.city_shapes:
                    raise ValueError(f"内部错误: 在 city_shapes 中找不到键 '{params['region_name']}'")
                city_info = self.city_shapes[params['region_name']]
                params['city_shape'] = city_info['shape']
                params['province_shape'] = self.province_shapes.get(city_info['province'])

                if not params['city_shape'] or not hasattr(params['city_shape'], 'bounds'):
                    raise ValueError(f"区域 '{params['region_name']}' 的地理形状无效或缺失！")

                west, south, east, north = params['city_shape'].bounds
                if not (north > south and east > west):
                    if params['province_shape'] and hasattr(params['province_shape'], 'bounds'):
                        p_west, p_south, p_east, p_north = params['province_shape'].bounds
                        if p_north > p_south and p_east > p_west:
                            west, south, east, north = p_west, p_south, p_east, p_north
                        else:
                            raise ValueError(f"区域 '{params['region_name']}' 及其关联省份的边界框均无效")
                    else:
                        raise ValueError(
                            f"区域 '{params['region_name']}' 的边界框无效: N{north}, S{south}, E{east}, W{west}")

                lon_buffer = (east - west) * 0.5
                lat_buffer = (north - south) * 0.5
                params['area_bounds'] = [west - lon_buffer, east + lon_buffer, south - lat_buffer, north + lat_buffer]

                lat_slice = slice(params['area_bounds'][2], params['area_bounds'][3])
                lon_slice = slice(params['area_bounds'][0], params['area_bounds'][1])
                ds_rect = self.ds.sel(latitude=lat_slice, longitude=lon_slice)

                if ds_rect.sizes['latitude'] == 0 or ds_rect.sizes['longitude'] == 0:
                    raise ValueError(
                        f"区域 '{params['region_name']}' 的扩展边界 [{lon_slice.start:.1f}-{lon_slice.stop:.1f}, {lat_slice.start:.1f}-{lat_slice.stop:.1f}] 完全落在数据范围之外。")

                mask = regionmask.mask_3D_geopandas(gpd.GeoSeries([params['city_shape']]), ds_rect.longitude,
                                                    ds_rect.latitude)
                ds_processed_masked = ds_rect.where(mask)

                first_var = self.data_vars_raw[0] if self.data_vars_raw else list(ds_processed_masked.data_vars)[0]
                if not self.data_vars_raw or ds_processed_masked[first_var].count().item() == 0:
                    self.root.after(0, self._update_status_safe,
                                    f"警告：区域 '{params['region_name']}' 精确范围内无数据点，已自动扩展为外接矩形分析。")
                    ds_processed = ds_rect
                else:
                    ds_processed = ds_processed_masked.isel(region=0)

                params['base_title_part'] = params['region_name'].replace('/', ' ')

            if not self.data_vars_raw or ds_processed[self.data_vars_raw[0]].count().item() == 0:
                raise ValueError(
                    f"空间区域 '{params['region_name']}' 内没有任何有效的气象数据点。\n\n这很可能是因为数据的空间分辨率过低，或选择的区域完全落在数据范围之外。")

            year = params['year']
            month = params['month']
            day = params['day']
            hour = params['hour']

            time_title_part = ""
            if params['chart_type'] in ["时间序列箱形图", "日历热力图", "统计直方图", "散布图"]:
                if not year: raise ValueError(f"绘制 {params['chart_type']} 需要在界面上明确选择一个年份。")
                ds_final = ds_processed.sel(time=ds_processed.time.dt.year == int(year))
                time_title_part = f"{year}年"
            else:
                time_mask = True
                if year: time_mask = time_mask & (
                        ds_processed.time.dt.year == int(year)); time_title_part += f"{year}年"
                if month: time_mask = time_mask & (
                        ds_processed.time.dt.month == int(month)); time_title_part += f"{month}月"
                if day: time_mask = time_mask & (ds_processed.time.dt.day == int(day)); time_title_part += f"{day}日"
                if hour: time_mask = time_mask & (
                        ds_processed.time.dt.hour == int(hour)); time_title_part += f"{hour}时"
                ds_final = ds_processed.sel(time=time_mask)
            params['time_title_part'] = time_title_part

            if ds_final.time.size == 0: raise ValueError("时间筛选后无数据点。")
            params['ds_final'] = ds_final

            var1 = params['var1']
            if not var1 or var1 not in ds_final: raise ValueError(
                f"选择的变量 '{params['var1_display']}' 无法在数据中找到。")

            if np.count_nonzero(~np.isnan(ds_final[params['var1']].values)) == 0:
                raise ValueError(f"数据有效性检查失败！在所选时空范围内，变量 '{params['var1']}' 的数值全部为空(NaN)。")

            self.root.after(0, self._execute_plotting, params)
        except Exception as e:
            traceback.print_exc()
            self.root.after(0, lambda exc=e: self._update_status_safe(f"数据准备失败：{exc}"))
            self.root.after(0, lambda exc=e: messagebox.showerror("错误", str(exc)))
            self.root.after(0, lambda: self.generate_button.config(state="normal"))

    def _execute_plotting(self, params):
        fig = None
        output_filename = None
        try:
            plt.close('all')
            self.root.after(0, self._update_status_safe, "数据准备完毕，正在主线程中生成图表...")

            if params['custom_title']:
                final_title = params['custom_title']
            else:
                var1_display_name = f"{VARIABLE_TRANSLATIONS.get(params['var1'], params['var1'])} ({params['var1']})"
                chart_type_cn = CHART_TYPE_TRANSLATIONS.get(params['chart_type'], "")
                base_title = f"{params['base_title_part']} {params['time_title_part']}"

                if params['chart_type'] == "日历热力图":
                    chart_type_cn = f"逐日{chart_type_cn}"

                final_title = f"{base_title}\n{var1_display_name} {chart_type_cn}"
                if params['chart_type'] == "风场图":
                    final_title = f"{base_title}\n{chart_type_cn}"
                elif params['chart_type'] == "散布图":
                    var2_display_name = f"{VARIABLE_TRANSLATIONS.get(params['var2'], params['var2'])} ({params['var2']})"
                    final_title = f"{base_title}\n{var1_display_name} 与 {var2_display_name} {chart_type_cn}"

            if params['custom_filename']:
                filename_base = params['custom_filename']
            else:
                timestamp = datetime.now().strftime("%Y%m%d-%H%M%S%f")[:-3]
                sanitized_title = re.sub(r'[\n/\\:\*\?\"<>\|]', '_', final_title.replace('\n', '_'))
                filename_base = f"{sanitized_title}_{timestamp}"
            output_filename = os.path.join(params['save_path'], f"{filename_base}.png")

            print(f"DEBUG: Executing plot for chart_type='{params['chart_type']}'")
            if params['chart_type'] == "时间段平均图 (地图)":
                fig = plot_geographical_map(params['ds_final'], params['var1'], final_title, params['area_bounds'],
                                            params['city_shape'], params['province_shape'], params['show_values'],
                                            params['dataset_name'])
            elif params['chart_type'] == "风场图":
                fig = plot_wind_map(params['ds_final'], params['var1'], params['var2'], final_title,
                                    params['area_bounds'], params['city_shape'], params['province_shape'],
                                    params['show_values'], params['dataset_name'])
            elif params['chart_type'] == "时间序列趋势图 (折线)":
                fig = plot_line_chart(params['ds_final'], params['var1'], final_title, params['aggregation'],
                                      params['dataset_name'])
            elif params['chart_type'] == "时间序列箱形图":
                fig = plot_boxplot(params['ds_final'], params['var1'], final_title, params['dataset_name'])
            elif params['chart_type'] == "日历热力图":
                fig = plot_calendar_heatmap(params['ds_final'], params['var1'], final_title, params['dataset_name'])
            elif params['chart_type'] == "统计直方图":
                fig = plot_histogram(params['ds_final'], params['var1'], final_title, params['dataset_name'])
            elif params['chart_type'] == "散布图":
                fig = plot_scatter(params['ds_final'], params['var1'], params['var2'], final_title,
                                   params['dataset_name'])
            else:
                raise ValueError(f"未知的图表类型或类型不匹配: '{params['chart_type']}'")

            if self.canvas_widget: self.canvas_widget.get_tk_widget().destroy()
            if self.toolbar_widget: self.toolbar_widget.destroy()
            if hasattr(self,
                       'placeholder_label') and self.placeholder_label.winfo_exists(): self.placeholder_label.destroy()
            self.canvas_widget = FigureCanvasTkAgg(fig, master=self.preview_frame)
            self.canvas_widget.draw()
            self.toolbar_widget = NavigationToolbar2Tk(self.canvas_widget, self.preview_frame)
            self.toolbar_widget.update()

            self.toolbar_widget.pack(side=tk.TOP, fill=tk.X)
            self.canvas_widget.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

            if params['save_mode'] == 'auto':
                fig.savefig(output_filename, dpi=300, bbox_inches='tight')
                plt.close(fig)
                self.root.after(50, self._update_status_safe, f"任务成功！图表已自动保存为：\n{output_filename}")
                self.save_chart_button.config(state="disabled")
                self.current_fig = None
                self.current_output_filename = None
            else:
                self.current_fig = fig
                self.current_output_filename = output_filename
                self.save_chart_button.config(state="normal")
                self.root.after(0, self._update_status_safe, f"图表已生成预览，请点击“保存当前图表”按钮进行保存。")

        except Exception as e:
            traceback.print_exc()
            self.root.after(0, lambda exc=e: self._update_status_safe(f"绘图失败：{exc}"))
            self.root.after(0, lambda exc=e: messagebox.showerror("绘图错误", str(exc)))
            if fig: plt.close(fig)
        finally:
            self.root.after(0, lambda: self.generate_button.config(state="normal"))

    def _save_current_chart(self):
        if self.current_fig and self.current_output_filename:
            try:
                self.current_fig.savefig(self.current_output_filename, dpi=300, bbox_inches='tight')
                plt.close(self.current_fig)
                self.root.after(50, self._update_status_safe,
                                f"图表已手动保存为：\n{self.current_output_filename}")
                self.current_fig = None
                self.current_output_filename = None
                self.save_chart_button.config(state="disabled")
            except Exception as e:
                traceback.print_exc()
                messagebox.showerror("保存失败", f"无法保存图表：\n{e}")
        else:
            messagebox.showwarning("无图表可保存", "请先生成一个图表，然后再进行保存。")
            self.save_chart_button.config(state="disabled")


# --- class LocalWeatherPlatform 结束 ---


if __name__ == "__main__":
    if gpd is None:
        messagebox.showerror("缺少关键工具库",
                             "未找到 'geopandas'。\n请在您的 Conda 环境中执行:\nconda install -c conda-forge geopandas regionmask")
        exit()
    try:
        import regionmask
    except ImportError:
        messagebox.showerror("缺少关键工具库",
                             "未找到 'regionmask'。\n请在您的 Conda 环境中执行:\nconda install -c conda-forge geopandas regionmask")
        exit()
    if OpenCC is None:
        messagebox.showerror("缺少可选库",
                             "未找到 'opencc-python-reimplemented' 库。\n繁体字地区名称将无法自动转换为简体。\n请运行 'pip install opencc-python-reimplemented' 安装。")
        # 不退出，允许程序继续运行，只是不做简繁转换

    root = tk.Tk()
    app = LocalWeatherPlatform(root)
    root.mainloop()