"""
多层薄膜光谱计算程序 - 图形用户界面版本
支持材料库选择、任意层数薄膜、多角度反射率计算

版本: 1.0
版权所有 © 2024 Jeffery & Cursor
保留所有权利
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext, simpledialog
import tkinter.font as tkfont
from material_database import MATERIALS, get_material_names, get_material_refractive_index, get_material_description, get_material_info,load_materials_from_file
import matplotlib
from manager_materials import MaterialManager
import json
import os
# 设置matplotlib使用中文字体
matplotlib.rcParams['font.sans-serif'] = ['SimHei',
                                          'Microsoft YaHei', 'SimSun', 'Arial Unicode MS']
matplotlib.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题


def calculate_multilayer_spectrum(
        n_layers,
        d_layers,
        wavelength_range=(
            400,
            800),
    incident_angle=0,
        polarization='s'):
    """计算多层薄膜的反射率光谱

    参数:
        n_layers: 各层折射率列表，包括入射介质和出射介质
        d_layers: 各层厚度列表（单位：nm），入射介质和出射介质厚度为0
        wavelength_range: 波长范围（单位：nm）
        incident_angle: 入射角度（单位：度）
        polarization: 偏振方式，'s'或'p'

    返回:
    wavelengths: 波长数组
    R: 反射率数组
    T: 透射率数组
    A: 吸收率数组
    """
    # 打印输入参数，帮助调试
    print(
        f"计算函数输入: 折射率={n_layers}, 厚度={d_layers}, 角度={incident_angle}°, 偏振={polarization}")

    # 生成波长数组
    wavelengths = np.linspace(wavelength_range[0], wavelength_range[1], 1000)

    # 初始化结果数组
    R = np.zeros_like(wavelengths, dtype=float)
    T = np.zeros_like(wavelengths, dtype=float)
    A = np.zeros_like(wavelengths, dtype=float)

    # 将入射角度转换为弧度
    theta_0 = np.deg2rad(incident_angle)

    # 特殊情况：只有基材没有膜层
    if len(n_layers) == 3 and len(d_layers) == 3:  # 入射介质、基材、出射介质
        # 使用菲涅耳公式直接计算反射率
        n1 = n_layers[0]  # 入射介质折射率
        n2 = n_layers[1]  # 基材折射率
        n3 = n_layers[2]  # 出射介质折射率

        # 计算入射角
        theta1 = theta_0

        # 计算基材中的角度（斯涅尔定律）
        try:
            theta2 = np.arcsin((n1 / n2) * np.sin(theta1))
        except BaseException:
            theta2 = 0

        # 计算出射角
        try:
            theta3 = np.arcsin((n2 / n3) * np.sin(theta2))
        except BaseException:
            theta3 = 0

        # 计算基材厚度（转换为nm）
        d = d_layers[1]

        # 对每个波长计算
        for i, wavelength in enumerate(wavelengths):
            # 计算空气-玻璃界面的反射率
            if polarization == 's':
                # s偏振的菲涅耳系数
                r12 = ((n1 * np.cos(theta1) - n2 * np.cos(theta2)) /
                       (n1 * np.cos(theta1) + n2 * np.cos(theta2)))
                r23 = ((n2 * np.cos(theta2) - n3 * np.cos(theta3)) /
                       (n2 * np.cos(theta2) + n3 * np.cos(theta3)))
            else:
                # p偏振的菲涅耳系数
                r12 = ((n2 * np.cos(theta1) - n1 * np.cos(theta2)) /
                       (n2 * np.cos(theta1) + n1 * np.cos(theta2)))
                r23 = ((n3 * np.cos(theta2) - n2 * np.cos(theta3)) /
                       (n3 * np.cos(theta2) + n2 * np.cos(theta3)))

            # 计算基材中的相位差
            delta = 2 * np.pi * n2 * d * np.cos(theta2) / wavelength

            # 考虑多重反射的干涉效应
            numerator = r12 + r23 * np.exp(2j * delta)
            denominator = 1 + r12 * r23 * np.exp(2j * delta)
            r = numerator / denominator

            # 计算反射率
            R[i] = np.abs(r) ** 2

            # 考虑基材的吸收
            if np.imag(n2) > 0:
                # 计算吸收系数 k = 4π·n_imag/λ
                k = 4 * np.pi * np.imag(n2) / wavelength
                # 计算透过率 T = exp(-k·d)
                absorption = 1 - np.exp(-k * d)
                A[i] = absorption * (1 - R[i])

            # 计算透射率（考虑能量守恒）
            T[i] = 1 - R[i] - A[i]

        print(f"使用菲涅耳公式计算单层基材的反射率")
        print(
            f"计算结果: 平均反射率={np.mean(R):.4f}, 最小={np.min(R):.4f}, 最大={np.max(R):.4f}")
        return wavelengths, R, T, A

    # 对每个波长计算反射率 - 多层膜情况
    for i, wavelength in enumerate(wavelengths):
        # 初始化吸收率
        A[i] = 0.0

        # 计算每层的角度
        thetas = [theta_0]  # 各层中的角度
        n0 = n_layers[0]  # 入射介质折射率

        for j in range(1, len(n_layers)):
            nj = n_layers[j]
            # 使用斯涅尔定律计算角度
            try:
                theta_j = np.arcsin((n0 / nj) * np.sin(theta_0))
                thetas.append(theta_j)
            except BaseException:
                # 处理全反射情况
                thetas.append(0)

        # 使用简单的传输矩阵法计算
        # 初始化传输矩阵
        M = np.array([[1, 0], [0, 1]], dtype=complex)

        # 从第一层到最后一层
        for j in range(1, len(n_layers)):
            # 计算当前层的折射率和厚度
            n_curr = n_layers[j]
            d_curr = d_layers[j]

            # 计算前一层的折射率
            n_prev = n_layers[j - 1]

            # 计算入射角和折射角
            theta_prev = thetas[j - 1]
            theta_curr = thetas[j]

            # 计算界面反射系数
            if polarization == 's':
                # s偏振的反射系数
                r = ((n_prev * np.cos(theta_prev) - n_curr * np.cos(theta_curr)) /
                     (n_prev * np.cos(theta_prev) + n_curr * np.cos(theta_curr)))

                # s偏振的透射系数
                t = (2 * n_prev * np.cos(theta_prev) / (n_prev *
                     np.cos(theta_prev) + n_curr * np.cos(theta_curr)))
            else:
                # p偏振的反射系数
                r = ((n_curr * np.cos(theta_prev) - n_prev * np.cos(theta_curr)) /
                     (n_curr * np.cos(theta_prev) + n_prev * np.cos(theta_curr)))

                # p偏振的透射系数
                t = (2 * n_prev * np.cos(theta_prev) / (n_curr *
                     np.cos(theta_prev) + n_prev * np.cos(theta_curr)))

            # 计算界面矩阵
            I = (1 / t) * np.array([[1, r], [r, 1]], dtype=complex)

            # 如果不是最后一层，还需要计算传播矩阵
            if j < len(n_layers) - 1:
                # 计算相位厚度
                beta = (
                    2 *
                    np.pi *
                    n_curr *
                    d_curr *
                    np.cos(theta_curr) /
                    wavelength)

                # 计算传播矩阵
                P = np.array([
                    [np.exp(-1j * beta), 0],
                    [0, np.exp(1j * beta)]
                ], dtype=complex)

                # 更新总传输矩阵
                M = M @ I @ P
            else:
                # 最后一层只需要界面矩阵
                M = M @ I

        # 计算反射系数
        r = M[1, 0] / M[0, 0]

        # 计算反射率
        R[i] = np.abs(r) ** 2

        # 添加安全检查，确保反射率在0到1之间
        if R[i] > 1:
            print(f"警告：波长{wavelength}nm处计算的反射率为{R[i]}，超出了物理范围，已修正为1")
            R[i] = 1.0
        elif R[i] < 0:
            print(f"警告：波长{wavelength}nm处计算的反射率为{R[i]}，低于物理下限，已修正为0")
            R[i] = 0.0

        # 计算透射率（考虑能量守恒）
        if np.imag(n_layers[-1]) == 0:  # 如果出射介质是无损耗的
            T[i] = 1 - R[i]
        else:
            T[i] = 0

        # 考虑吸收
        for j in range(1, len(n_layers) - 1):
            if np.imag(n_layers[j]) > 0:
                # 计算吸收系数
                k = 4 * np.pi * np.imag(n_layers[j]) / wavelength
                # 计算透过率
                absorption = 1 - np.exp(-k * d_layers[j])
                # 累积到总吸收率
                A[i] += absorption * T[i]

        # 重新计算透射率（考虑吸收）
        T[i] = 1 - R[i] - A[i]

    # 打印计算结果的统计信息
    print(
        f"计算结果: 平均反射率={np.mean(R):.4f}, 最小={np.min(R):.4f}, 最大={np.max(R):.4f}")

    return wavelengths, R, T, A


class LayerFrame(ttk.Frame):
    """表示单个薄膜层的框架"""

    def __init__(
            self,
            parent,
            layer_num,
            is_medium=False,
            is_substrate=False,
            is_exit_medium=False,
            *args,
            **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.parent = parent
        self.layer_num = layer_num
        self.is_medium = is_medium
        self.is_substrate = is_substrate
        self.is_exit_medium = is_exit_medium

        # 获取材料列表
        self.material_names = get_material_names()

        # 添加标签属性
        self.label = None

        # 创建控件
        self.create_widgets()

    def create_widgets(self):
        # 层标签
        if self.is_medium:
            label_text = "入射介质:"
        elif self.is_exit_medium:
            label_text = "出射介质:"
        elif self.is_substrate:
            label_text = "基材材料:"
        else:
            label_text = f"第 {self.layer_num} 层:"

        self.label = ttk.Label(self, text=label_text, width=10)
        self.label.grid(row=0, column=0, padx=2, pady=2, sticky="w")

        # 材料选择下拉框
        self.material_var = tk.StringVar()
        if self.is_medium:
            default_material = "Air"
        elif self.is_exit_medium:
            default_material = "Air"
        elif self.is_substrate:
            default_material = "Glass"
        else:
            default_material = "SiO2"
        self.material_var.set(default_material)

        self.material_combo = ttk.Combobox(
            self,
            textvariable=self.material_var,
            values=self.material_names,
            width=12)
        self.material_combo.grid(row=0, column=1, padx=2, pady=2, sticky="w")
        self.material_combo.bind(
            "<<ComboboxSelected>>",
            self.on_material_selected)

        # 折射率输入框
        ttk.Label(
            self,
            text="折射率:",
            width=6).grid(
            row=0,
            column=2,
            padx=2,
            pady=2,
            sticky="w")
        self.n_var = tk.StringVar()
        n_value = get_material_refractive_index(default_material)
        self.n_var.set(str(n_value))
        self.n_entry = ttk.Entry(self, textvariable=self.n_var, width=18)
        self.n_entry.grid(row=0, column=3, padx=2, pady=2, sticky="w")

        # 厚度输入框
        if not self.is_medium and not self.is_exit_medium:
            ttk.Label(
                self,
                text="厚度(nm):",
                width=8).grid(
                row=0,
                column=4,
                padx=2,
                pady=2,
                sticky="w")
            self.d_var = tk.StringVar()
            self.d_var.set("100")  # 默认厚度
            self.d_entry = ttk.Entry(self, textvariable=self.d_var, width=10)
            self.d_entry.grid(row=0, column=5, padx=2, pady=2, sticky="w")
        else:
            self.d_var = tk.StringVar()
            self.d_var.set("0")  # 入射介质和出射介质厚度为0

    def on_material_selected(self, event):
        """当材料选择改变时更新折射率"""
        material = self.material_var.get()
        n_value = get_material_refractive_index(material)

        if material == "Custom":
            self.n_entry.config(state="normal")
            self.n_var.set("")
        else:
            self.n_entry.config(state="normal")
            if n_value is not None:
                self.n_var.set(str(n_value))

    def get_values(self):
        """获取当前层的材料、折射率和厚度"""
        material = self.material_var.get()

        # 解析折射率（可能是复数）
        n_str = self.n_var.get()
        try:
            if '+' in n_str or '-' in n_str and 'j' in n_str:
                n = complex(n_str)
            else:
                n = float(n_str)
        except ValueError:
            messagebox.showerror("输入错误", f"无效的折射率: {n_str}")
            return None, None, None

        # 解析厚度
        d_str = self.d_var.get()
        try:
            d = float(d_str)
        except ValueError:
            messagebox.showerror("输入错误", f"无效的厚度: {d_str}")
            return None, None, None

        return material, n, d

    def update_label(self, layer_num):
        """更新层号标签"""
        self.label.configure(text=f"第 {layer_num} 层:")


class PlotWindow(tk.Toplevel):
    """独立的图表窗口"""

    def __init__(self, parent, title="多层薄膜光谱"):
        super().__init__(parent)
        self.title(title + " v1.0")
        self.geometry("1000x700")

        # 创建图表
        self.fig = Figure(figsize=(10, 7), dpi=100)
        self.ax = self.fig.add_subplot(111)

        # 将图形嵌入窗口
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)

        # 添加工具栏
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()
        
        # 创建固定在左上角的文本框用于显示数据
        self.info_text = self.fig.text(0.02, 0.98, '', 
                                      transform=self.fig.transFigure,
                                      verticalalignment='top',
                                      bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7),
                                      fontsize=10)
        self.info_text.set_visible(False)
        
        # 添加垂直线用于标记当前鼠标位置
        self.vline = None
        
        # 存储曲线数据
        self.line_data = []
        
        # 使用Tkinter的事件绑定
        self.canvas_widget.bind("<Motion>", self.on_mouse_move_tk)
        
        # 设置关闭窗口的行为
        self.protocol("WM_DELETE_WINDOW", self.on_close)
        
    def on_mouse_move_tk(self, event):
        """处理Tkinter的鼠标移动事件"""
        # 将Tkinter坐标转换为Matplotlib坐标
        x, y = event.x, event.y
        
        # 转换为数据坐标
        try:
            # 获取图表区域的边界
            bbox = self.ax.get_position()
            canvas_width = self.canvas_widget.winfo_width()
            canvas_height = self.canvas_widget.winfo_height()
            
            # 检查鼠标是否在图表区域内
            ax_x0 = bbox.x0 * canvas_width
            ax_y0 = bbox.y0 * canvas_height
            ax_x1 = bbox.x1 * canvas_width
            ax_y1 = bbox.y1 * canvas_height
            
            if not (ax_x0 <= x <= ax_x1 and ax_y0 <= y <= ax_y1):
                # 鼠标不在图表区域内
                self.info_text.set_visible(False)
                if self.vline:
                    self.vline.set_visible(False)
                self.canvas.draw_idle()
                return
            
            # 计算相对位置
            rel_x = (x - ax_x0) / (ax_x1 - ax_x0)
            rel_y = 1.0 - (y - ax_y0) / (ax_y1 - ax_y0)  # 反转Y轴
            
            # 转换为数据坐标
            x_min, x_max = self.ax.get_xlim()
            y_min, y_max = self.ax.get_ylim()
            
            data_x = x_min + rel_x * (x_max - x_min)
            data_y = y_min + rel_y * (y_max - y_min)
            
            # 更新或创建垂直线
            if self.vline is None:
                self.vline = self.ax.axvline(data_x, color='gray', linestyle='--', alpha=0.5)
            else:
                # 正确设置垂直线位置 - 使用列表
                self.vline.set_xdata([data_x, data_x])
                self.vline.set_visible(True)
            
            # 显示所有曲线在当前波长的反射率
            self.show_all_reflectance_at_wavelength(data_x)
            
        except Exception as e:
            print(f"坐标转换错误: {e}")
            self.info_text.set_visible(False)
            if self.vline:
                self.vline.set_visible(False)
            self.canvas.draw_idle()

    def show_all_reflectance_at_wavelength(self, x_mouse):
        """显示所有曲线在指定波长处的反射率值"""
        if not self.line_data:
            return
            
        info_text = f"波长: {x_mouse:.1f} nm\n"
        
        # 对于每条曲线，找到最接近的波长点
        for wavelengths, reflectance, label in self.line_data:
            # 找到最接近的波长索引
            idx = min(range(len(wavelengths)), key=lambda i: abs(wavelengths[i] - x_mouse))
            wavelength = wavelengths[idx]
            r_value = reflectance[idx] * 100  # 转换为百分比
            
            # 添加到显示文本
            info_text += f"{label}: {r_value:.2f}%\n"
        
        # 更新文本框内容
        self.info_text.set_text(info_text)
        self.info_text.set_visible(True)
        
        # 重绘画布
        self.canvas.draw_idle()

    def on_close(self):
        """关闭窗口时的行为"""
        self.withdraw()  # 隐藏窗口而不是销毁它

    def clear_plot(self):
        """清除当前图表"""
        self.ax.clear()
        self.ax.set_xlabel('波长 (nm)')
        self.ax.set_ylabel('反射率 (%)')
        self.ax.set_title('反射率光谱')
        self.ax.grid(True)
        self.fig.tight_layout()
        
        # 重置垂直线
        self.vline = None
        
        # 清空曲线数据
        self.line_data = []
        
        # 隐藏信息文本
        self.info_text.set_visible(False)
        
        self.canvas.draw()

    def plot_spectrum(self, wavelengths, reflectance, polarization, angle):
        """绘制一条反射率曲线"""
        # 定义颜色和线型
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        linestyles = {'s': '-', 'p': '--'}

        color = colors[int(angle / 15) % len(colors)]
        linestyle = linestyles[polarization]
        # 将反射率转换为百分比
        reflectance_percent = reflectance * 100
        self.ax.plot(
            wavelengths,
            reflectance_percent,
            color=color,
            linestyle=linestyle,
            linewidth=2,
            label=f'{polarization}偏振 {angle}°')
            
        # 存储曲线数据用于鼠标交互
        self.line_data.append((wavelengths, reflectance, f'{polarization}偏振 {angle}°'))

        self.ax.set_xlabel('波长 (nm)')
        self.ax.set_ylabel('反射率 (%)')
        self.ax.set_title('反射率光谱')
        self.ax.grid(True)

        # 设置x轴范围
        self.ax.set_xlim(min(wavelengths), max(wavelengths))

        # 设置y轴范围
        max_R = max(reflectance_percent)
        self.ax.set_ylim(0, max(105, max_R * 1.1))

        # 添加图例
        self.ax.legend()

        # 更新图表
        self.fig.tight_layout()
        self.canvas.draw()

        # 显示窗口
        self.deiconify()


class MultilayerFilmApp(tk.Tk):
    """多层薄膜光谱计算程序的主应用"""

    def __init__(self):
        super().__init__()

        self.material_manager = MaterialManager()
        load_materials_from_file("custom_materials.json")
        self.material_manager = MaterialManager(self)

        # 设置窗口标题和大小
        self.title("多层薄膜光谱计算程序")
        self.geometry("1200x800")

        # 创建变量
        self.incident_medium_var = tk.StringVar(value="Air")
        self.substrate_var = tk.StringVar(value="Glass")
        self.substrate_thickness_var = tk.StringVar(value="1")
        self.substrate_unit_var = tk.StringVar(value="mm")
        self.exit_medium_var = tk.StringVar(value="Air")

        # 计算参数变量
        self.wavelength_var = tk.StringVar(value="400-800")
        self.angle_var = tk.StringVar(value="0")  # 默认入射角度为0度
        self.polarization_s_var = tk.BooleanVar(value=True)   # 默认选择s光
        self.polarization_p_var = tk.BooleanVar(value=False)  # 默认不选择p光

        # 初始化层列表
        self.layer_frames = []
        self.exit_layer_frames = []

        # 创建主框架，修改pack参数，删除expand=True
        self.main_frame = ttk.Frame(self)
        self.main_frame.pack(fill="x", padx=10, pady=(10, 0))  # 只在顶部添加padding

        # 创建顶部框架（用于介质选择）
        self.top_frame = ttk.Frame(self.main_frame)
        self.top_frame.pack(fill="x", pady=(0, 5))  # 减小顶部padding

        # 创建中部框架（用于膜层设置）
        self.middle_frame = ttk.Frame(self.main_frame)
        self.middle_frame.pack(fill="x", pady=5)  # 改为fill="x"

        # 创建底部框架（用于控制按钮和设置）
        self.bottom_frame = ttk.Frame(self.main_frame)
        self.bottom_frame.pack(fill="x", pady=5)

        # 创建控件
        self.create_input_widgets()

        # 创建图表窗口（初始隐藏）
        self.plot_window = PlotWindow(self)
        self.plot_window.withdraw()

        # 默认添加标准的双面增透膜结构
        self.apply_preset("双面增透膜")

    def create_input_widgets(self):
        """创建输入区域的控件"""
        # 顶部控制区域 - 使用垂直布局
        top_control = ttk.Frame(self.main_frame)
        top_control.pack(fill="x", pady=5)

        # 修改顶部按钮布局 - 第一行包含左边材料管理、中间计算结果输出框、右边计算光谱
        button_frame = ttk.Frame(top_control)
        button_frame.pack(fill="x", padx=10, pady=5)

        # 左侧放材料管理按钮
        ttk.Button(button_frame, text="材料管理", command=self.manage_materials,
                   style='Big.TButton').pack(side="left", padx=(0, 5))

        # 右侧放计算光谱按钮
        ttk.Button(button_frame, text="计算光谱", command=self.calculate_spectrum,
                   style='Big.TButton').pack(side="right", padx=(5, 0))

        # 中间放计算结果输出框
        output_frame = ttk.LabelFrame(button_frame, text="计算结果")
        output_frame.pack(side="left", fill="both", expand=True, padx=5)

        # 创建文本显示区域
        self.output_text = scrolledtext.ScrolledText(
            output_frame, height=6, width=50)
        self.output_text.pack(fill="both", expand=True, padx=5, pady=5)

        # 第二行：预设膜层结构
        preset_frame = ttk.LabelFrame(top_control, text="预设膜层结构")
        preset_frame.pack(fill="x", pady=5)

        preset_content = ttk.Frame(preset_frame)
        preset_content.pack(fill="x", padx=10, pady=5)

        presets = [("双面增透膜", "double_side_ar"), ("四分之一波长高反射膜", "quarter_wave_stack"),
                   ("宽带增透膜", "broadband_ar"), ("单层增透膜", "single_layer_ar"),
                   ("金属反射膜", "metal_reflector"),
                   ("低通高反膜", "low_pass_filter")]

        ttk.Label(preset_content, text="选择预设结构:").pack(side="left", padx=5)
        self.preset_var = tk.StringVar(value="四分之一波长高反射膜")
        preset_combo = ttk.Combobox(
            preset_content, textvariable=self.preset_var, values=[
                p[0] for p in presets], width=25)
        preset_combo.pack(side="left", padx=5)
        ttk.Button(
            preset_content,
            text="应用预设",
            command=lambda: self.apply_preset(
                self.preset_var.get())).pack(
            side="left",
            padx=5)

        # 第三行：将三个介质放在同一行
        media_frame = ttk.LabelFrame(top_control, text="介质设置")
        media_frame.pack(fill="x", pady=5)

        media_content = ttk.Frame(media_frame)
        media_content.pack(fill="x", padx=10, pady=5)

        # 入射介质
        incident_frame = ttk.Frame(media_content)
        incident_frame.pack(side="left", padx=20)
        ttk.Label(incident_frame, text="入射介质:").pack(side="left", padx=2)
        self.incident_medium_combo = ttk.Combobox(
            incident_frame,
            textvariable=self.incident_medium_var,
            values=get_material_names(),
            width=15)
        self.incident_medium_combo.pack(side="left", padx=2)

        # 基材设置
        substrate_frame = ttk.Frame(media_content)
        substrate_frame.pack(side="left", padx=20)
        ttk.Label(substrate_frame, text="基材:").pack(side="left", padx=2)
        self.substrate_combo = ttk.Combobox(
            substrate_frame,
            textvariable=self.substrate_var,
            values=get_material_names(),
            width=15)
        self.substrate_combo.pack(side="left", padx=2)
        ttk.Label(substrate_frame, text="厚度:").pack(side="left", padx=2)
        self.substrate_thickness_entry = ttk.Entry(
            substrate_frame, textvariable=self.substrate_thickness_var, width=8)
        self.substrate_thickness_entry.pack(side="left", padx=2)
        self.substrate_unit_combo = ttk.Combobox(
            substrate_frame, textvariable=self.substrate_unit_var, values=[
                "nm", "μm", "mm"], width=5, state="readonly")
        self.substrate_unit_combo.pack(side="left", padx=2)

        # 出射介质
        exit_frame = ttk.Frame(media_content)
        exit_frame.pack(side="left", padx=20)
        ttk.Label(exit_frame, text="出射介质:").pack(side="left", padx=2)
        self.exit_medium_combo = ttk.Combobox(
            exit_frame,
            textvariable=self.exit_medium_var,
            values=get_material_names(),
            width=15)
        self.exit_medium_combo.pack(side="left", padx=2)

        # 第六行：计算参数
        calc_frame = ttk.LabelFrame(top_control, text="计算参数")
        calc_frame.pack(fill="x", pady=5)

        calc_content = ttk.Frame(calc_frame)
        calc_content.pack(fill="x", padx=10, pady=5)

        # 使用水平布局排列三个参数
        # 波长范围
        wavelength_frame = ttk.Frame(calc_content)
        wavelength_frame.pack(side="left", padx=20)
        ttk.Label(wavelength_frame, text="波长范围:").pack(side="left", padx=2)
        ttk.Entry(
            wavelength_frame,
            textvariable=self.wavelength_var,
            width=15).pack(
            side="left",
            padx=2)
        ttk.Label(wavelength_frame, text="nm").pack(side="left", padx=2)

        # 入射角度
        angle_frame = ttk.Frame(calc_content)
        angle_frame.pack(side="left", padx=20)
        ttk.Label(angle_frame, text="入射角度:").pack(side="left", padx=2)
        ttk.Entry(
            angle_frame,
            textvariable=self.angle_var,
            width=15).pack(
            side="left",
            padx=2)
        ttk.Label(angle_frame, text="°（多个角度用逗号分隔）").pack(side="left", padx=2)

        # 偏振
        polarization_frame = ttk.Frame(calc_content)
        polarization_frame.pack(side="left", padx=20)
        ttk.Label(polarization_frame, text="偏振:").pack(side="left", padx=2)
        ttk.Checkbutton(
            polarization_frame,
            text="s光",
            variable=self.polarization_s_var).pack(
            side="left",
            padx=2)
        ttk.Checkbutton(
            polarization_frame,
            text="p光",
            variable=self.polarization_p_var).pack(
            side="left",
            padx=2)

        # 下部分膜层设置区域（左右两列）
        films_frame = ttk.Frame(self.main_frame)
        films_frame.pack(fill="both", expand=True, pady=5)

        # 左侧上层膜层框架修改
        left_frame = ttk.LabelFrame(films_frame, text="上层膜层设置")
        left_frame.pack(side="left", fill="both", expand=True, padx=5)

        # 添加/删除层按钮 - 移到顶部
        buttons_frame = ttk.Frame(left_frame)
        buttons_frame.pack(fill="x", pady=2)
        ttk.Button(
            buttons_frame,
            text="添加膜层",
            command=self.add_layer).pack(
            side="left",
            padx=5)
        ttk.Button(
            buttons_frame,
            text="删除膜层",
            command=self.remove_layer).pack(
            side="left",
            padx=5)

        # 膜层显示区域
        self.layers_canvas = tk.Canvas(left_frame, height=400)
        self.layers_scrollbar = ttk.Scrollbar(
            left_frame, orient="vertical", command=self.layers_canvas.yview)

        self.layers_frame = ttk.Frame(self.layers_canvas)
        self.layers_frame.bind(
            "<Configure>", lambda e: self.layers_canvas.configure(
                scrollregion=self.layers_canvas.bbox("all")))

        # 基材标签放在最下面
        style = ttk.Style()
        style.configure('Substrate.TLabel', font=('SimHei', 10, 'bold'))
        self.substrate_label = ttk.Label(
            self.layers_frame, text="基材", style='Substrate.TLabel')
        self.substrate_label.pack(side="bottom", pady=5)

        self.layers_canvas.create_window(
            (0, 0), window=self.layers_frame, anchor="nw")
        self.layers_canvas.configure(yscrollcommand=self.layers_scrollbar.set)

        self.layers_scrollbar.pack(side="right", fill="y")
        self.layers_canvas.pack(side="left", fill="both", expand=True)

        # 右侧下层膜层
        right_frame = ttk.LabelFrame(films_frame, text="下层膜层设置")
        right_frame.pack(side="left", fill="both", expand=True, padx=5)

        # 添加/删除基材下方膜层按钮 - 移到顶部
        exit_buttons_frame = ttk.Frame(right_frame)
        exit_buttons_frame.pack(fill="x", pady=2)

        # 基材标签放在按钮下方
        style.configure('SubstrateDown.TLabel', font=('SimHei', 10, 'bold'))
        substrate_down_label = ttk.Label(
            right_frame, text="基材", style='SubstrateDown.TLabel')
        substrate_down_label.pack(pady=5)

        ttk.Button(
            exit_buttons_frame,
            text="添加膜层",
            command=self.add_exit_layer).pack(
            side="left",
            padx=5)
        ttk.Button(
            exit_buttons_frame,
            text="删除膜层",
            command=self.remove_exit_layer).pack(
            side="left",
            padx=5)

        # 膜层显示区域
        self.exit_layers_canvas = tk.Canvas(right_frame, height=400)
        self.exit_layers_scrollbar = ttk.Scrollbar(
            right_frame, orient="vertical", command=self.exit_layers_canvas.yview)

        self.exit_layers_frame = ttk.Frame(self.exit_layers_canvas)
        self.exit_layers_frame.bind(
            "<Configure>", lambda e: self.exit_layers_canvas.configure(
                scrollregion=self.exit_layers_canvas.bbox("all")))

        self.exit_layers_canvas.create_window(
            (0, 0), window=self.exit_layers_frame, anchor="nw")
        self.exit_layers_canvas.configure(
            yscrollcommand=self.exit_layers_scrollbar.set)

        self.exit_layers_scrollbar.pack(side="right", fill="y")
        self.exit_layers_canvas.pack(side="left", fill="both", expand=True)

    def add_layer(self):
        """添加一个新的膜层"""
        try:
            # 打印调试信息
            print("开始添加新膜层")

            # 获取现有膜层数量
            layer_num = len(self.layer_frames) + 1

            # 创建新的膜层框架
            layer_frame = LayerFrame(self.layers_frame, layer_num)

            # 将新膜层添加到列表开头，这样最新的膜层会显示在最上面
            self.layer_frames.insert(0, layer_frame)

            # 重新排列所有膜层，确保它们按照从基材向上的顺序显示
            for frame in self.layer_frames:
                frame.pack_forget()  # 先移除所有膜层

            # 从第一层开始重新添加，确保顺序正确
            for i, frame in enumerate(reversed(self.layer_frames)):
                frame.update_label(i + 1)  # 更新层号标签
                frame.pack(fill="x", pady=2, before=self.substrate_label)

            print(f"成功添加膜层，当前共有{len(self.layer_frames)}层")

            # 更新界面
            self.update()
            return True
        except Exception as e:
            print(f"添加膜层时发生错误: {str(e)}")
            import traceback
            traceback.print_exc()
            return False

    def remove_layer(self):
        """删除最上面的膜层"""
        try:
            if not self.layer_frames:
                print("没有可删除的膜层")
                return False

            # 获取最上面的膜层（列表中的第一个元素）并销毁
            layer_frame = self.layer_frames.pop(0)
            layer_frame.destroy()

            # 更新剩余膜层的标签
            for i, frame in enumerate(reversed(self.layer_frames)):
                frame.update_label(i + 1)

            print(f"成功删除膜层，当前剩余{len(self.layer_frames)}层")

            # 更新界面
            self.update()
            return True
        except Exception as e:
            print(f"删除膜层时发生错误: {str(e)}")
            import traceback
            traceback.print_exc()
            return False

    def add_exit_layer(self):
        """添加新的基材下方膜层"""
        exit_layer_frame = LayerFrame(
            self.exit_layers_frame, len(
                self.exit_layer_frames) + 1)  # 从1开始
        exit_layer_frame.pack(fill="x", pady=2)
        self.exit_layer_frames.append(exit_layer_frame)  # 添加到下方膜层列表

        # 更新滚动区域
        self.exit_layers_canvas.configure(
            scrollregion=self.exit_layers_canvas.bbox("all"))

    def remove_exit_layer(self):
        """删除基材下方的膜层"""
        if len(self.exit_layer_frames) == 0:  # 没有膜层可以删除
            messagebox.showinfo("提示", "没有更多的基材下方膜层可以删除")
            return

        # 移除最上面的基材下方膜层
        layer_to_remove = self.exit_layer_frames.pop()  # 移除最后一层
        layer_to_remove.destroy()

        # 更新滚动区域
        self.exit_layers_canvas.configure(
            scrollregion=self.exit_layers_canvas.bbox("all"))

    def on_substrate_selected(self, event):
        """当基材选择改变时更新折射率"""
        material = self.substrate_var.get()
        n_value = get_material_refractive_index(material)

        # 不再修改厚度值，只更新折射率显示
        if n_value is not None:
            # 这里可以添加折射率显示的代码，如果需要的话
            pass

    def calculate_spectrum(self):
        """计算并显示光谱"""
        try:
            # 收集所有层的参数
            n_layers = []
            d_layers = []
            
            # 获取入射介质折射率
            incident_material = self.incident_medium_var.get()
            incident_n = get_material_refractive_index(incident_material)
            if incident_n is None:
                messagebox.showerror("输入错误", f"无法获取入射介质 {incident_material} 的折射率")
                return
            n_layers.append(incident_n)
            d_layers.append(0)  # 入射介质厚度为0
            
            # 获取上方膜层参数 - 修改这里，反转顺序以匹配界面显示
            # 界面从上到下显示的是从入射介质到基材的顺序
            for layer_frame in reversed(self.layer_frames):  # 使用reversed确保顺序正确
                material, n, d = layer_frame.get_values()
                if material is None:  # 输入错误
                    return
                
                n_layers.append(n)
                d_layers.append(d)
            
            # 获取基材折射率
            substrate_material = self.substrate_var.get()
            substrate_n = get_material_refractive_index(substrate_material)
            if substrate_n is None:
                messagebox.showerror("输入错误", f"无法获取基材 {substrate_material} 的折射率")
                return
            n_layers.append(substrate_n)
            
            # 获取基材厚度并转换为纳米
            try:
                thickness_substrate = float(self.substrate_thickness_var.get())
                if thickness_substrate <= 0:
                    messagebox.showerror("输入错误", "基材厚度必须大于0")
                    return
                    
                if self.substrate_unit_var.get() == "μm":
                    thickness_substrate *= 1000  # 转换为nm
                elif self.substrate_unit_var.get() == "mm":
                    thickness_substrate *= 1000000  # 转换为nm
                    
                # 添加合理性检查
                if thickness_substrate < 100:  # 小于100nm可能不合理
                    if not messagebox.askyesno("警告", 
                        "基材厚度似乎过小（小于100nm），是否继续计算？"):
                        return
                
                d_layers.append(thickness_substrate)
            except ValueError:
                messagebox.showerror("输入错误", "基材厚度必须是有效的数字")
                return
            
            # 获取下方膜层参数
            for layer_frame in self.exit_layer_frames:
                material, n, d = layer_frame.get_values()
                if material is None:  # 输入错误
                    return
                
                n_layers.append(n)
                d_layers.append(d)
            
            # 获取出射介质折射率
            exit_material = self.exit_medium_var.get()
            exit_n = get_material_refractive_index(exit_material)
            if exit_n is None:
                messagebox.showerror("输入错误", f"无法获取出射介质 {exit_material} 的折射率")
                return
            n_layers.append(exit_n)
            d_layers.append(0)  # 出射介质厚度为0
            
            # 获取波长范围
            wavelength_range = self.wavelength_var.get().split('-')
            if len(wavelength_range) != 2:
                messagebox.showerror("输入错误", "波长范围格式应为'最小值-最大值'")
                return
            
            try:
                min_wavelength = float(wavelength_range[0])
                max_wavelength = float(wavelength_range[1])
                
                if min_wavelength >= max_wavelength:
                    messagebox.showerror("输入错误", "波长范围的最小值必须小于最大值")
                    return
                    
                if min_wavelength <= 0 or max_wavelength <= 0:
                    messagebox.showerror("输入错误", "波长必须大于0")
                    return
            except ValueError:
                messagebox.showerror("输入错误", "波长范围必须是有效的数字")
                return
            
            # 获取入射角度
            try:
                # 清理输入，移除空格
                angle_input = self.angle_var.get().strip()
                print(f"原始输入: {angle_input}")
                
                # 替换中文逗号为英文逗号
                angle_input = angle_input.replace('，', ',')
                print(f"清理后输入: {angle_input}")
                
                # 解析角度列表
                if ',' in angle_input:
                    # 多个角度，用逗号分隔
                    incident_angles = [float(angle.strip()) for angle in angle_input.split(',') if angle.strip()]
                else:
                    # 单个角度
                    incident_angles = [float(angle_input)] if angle_input else [0.0]
                
                print(f"解析后的角度列表: {incident_angles}")
                
                if not incident_angles:
                    messagebox.showerror("输入错误", "请输入至少一个有效的入射角度")
                    return
                
                # 检查角度是否在合理范围内
                for angle in incident_angles:
                    if angle < 0 or angle >= 90:
                        messagebox.showerror("输入错误", f"入射角度必须在0到90度之间，当前角度：{angle}")
                        return
                        
            except Exception as e:
                messagebox.showerror("输入错误", f"角度输入格式错误：{str(e)}\n请输入数字，多个角度用逗号分隔")
                return
            
            # 获取偏振状态
            polarizations = []
            if self.polarization_s_var.get():
                polarizations.append('s')
            if self.polarization_p_var.get():
                polarizations.append('p')
            
            if not polarizations:
                messagebox.showerror("输入错误", "请至少选择一种偏振状态")
                return
            
            # 打印调试信息
            print(f"计算参数: 层数={len(n_layers)}, 偏振={','.join(polarizations)}, 入射角={','.join(map(str, incident_angles))}")
            print(f"折射率列表: {n_layers}")
            print(f"厚度列表: {d_layers}")
            
            # 在绘制新图表前先清除之前的图表
            self.plot_window.clear_plot()
            
            # 清空输出文本框
            self.output_text.delete(1.0, tk.END)
            
            # 添加计算结果标题
            self.output_text.insert(tk.END, "计算结果摘要：\n", "title")
            self.output_text.tag_configure("title", font=("SimHei", 10, "bold"))
            
            # 计算并绘制光谱
            for angle in incident_angles:
                for polarization in polarizations:
                    wavelengths, R, _, _ = calculate_multilayer_spectrum(
                        n_layers, 
                        d_layers, 
                        wavelength_range=(min_wavelength, max_wavelength),
                        incident_angle=angle,
                        polarization=polarization
                    )
                    
                    # 绘制光谱
                    self.plot_window.plot_spectrum(wavelengths, R, polarization, angle)
                    
                    # 计算平均反射率、最小值和最大值
                    avg_R = np.mean(R) * 100  # 转换为百分比
                    min_R = np.min(R) * 100
                    max_R = np.max(R) * 100
                    
                    # 添加到输出文本框
                    result_text = f"{polarization}偏振, {angle}°角度: 平均反射率 = {avg_R:.2f}%, 范围 = [{min_R:.2f}% - {max_R:.2f}%]\n"
                    self.output_text.insert(tk.END, result_text)
                    
                    # 同时在终端中打印结果
                    print(f"结果统计: {polarization}偏振, {angle}°角度: 平均反射率 = {avg_R:.2f}%, 范围 = [{min_R:.2f}% - {max_R:.2f}%]")
                
            # 检查计算结果
            if len(wavelengths) == 0 or len(R) == 0:
                messagebox.showerror("计算错误", "计算结果为空")
                return
                
            print(f"计算完成: 波长点数={len(wavelengths)}, 反射率点数={len(R)}")
            print(f"反射率范围: {np.min(R)}-{np.max(R)}")
            
            # 显示图表窗口
            self.plot_window.deiconify()
            
        except Exception as e:
            print(f"计算过程中发生错误: {str(e)}")
            import traceback
            traceback.print_exc()
            messagebox.showerror("计算错误", f"计算过程中发生错误: {str(e)}")

    def apply_preset(self, preset_name):
        """应用预设膜层结构"""
        try:
            print(f"正在应用预设: {preset_name}")

            # 先清除现有膜层
            while self.layer_frames:
                success = self.remove_layer()
                if not success:
                    print("清除现有膜层失败")
                    break

            # 清空出射层帧
            for layer_frame in self.exit_layer_frames:
                layer_frame.destroy()
            self.exit_layer_frames.clear()

            # 添加调试信息
            print(f"已清除所有现有膜层，开始应用新预设: {preset_name}")

            # 设计波长 (nm) - 所有预设共用
            design_wavelength = 550

            # 选择预设 - 统一处理格式，避免不同名称无法匹配
            preset_name_lower = preset_name.lower()

            # 双面增透膜 (Double-side AR coating)
            if "双面增透膜" in preset_name or "double_side_ar" in preset_name_lower:
                print("正在应用双面增透膜预设...")

                # 设置入射介质
                self.incident_medium_var.set("Air")

                # 使用MgF2作为增透材料
                material = "MgF2"

                # 固定厚度值，而不是计算 - 使用典型值
                thickness = 100.0  # nm

                # 添加上表面膜层 - 从入射介质到基材方向
                self.add_layer()
                if self.layer_frames:
                    layer = self.layer_frames[0]
                    layer.material_var.set(material)
                    layer.d_var.set(f"{thickness:.1f}")
                    layer.on_material_selected(None)

                # 基材设为玻璃
                self.substrate_var.set("Glass")
                self.substrate_thickness_var.set("1")
                self.substrate_unit_var.set("mm")
                if hasattr(self, 'on_substrate_selected'):
                    self.on_substrate_selected(None)

                # 出射介质设为空气
                self.exit_medium_var.set("Air")

                # 添加下表面膜层 - 从基材到出射介质方向
                self.add_exit_layer()
                if self.exit_layer_frames:
                    layer = self.exit_layer_frames[0]
                    layer.material_var.set(material)
                    layer.d_var.set(f"{thickness:.1f}")
                    layer.on_material_selected(None)

            # 1. 四分之一波长高反射膜 (Quarter-wave stack)
            elif "四分之一波长高反射膜" in preset_name or "quarter_wave" in preset_name_lower:
                print("正在应用四分之一波长高反射膜预设...")

                # 设置入射介质
                self.incident_medium_var.set("Air")

                # 高低折射率材料
                high_index_material = "TiO2"
                low_index_material = "SiO2"

                # 直接使用固定厚度值，而不是计算
                high_thickness = 57.0  # nm
                low_thickness = 94.0  # nm

                # 创建7层高反射膜 (HLHLHLH) - 从入射介质到基材方向
                materials = [
                    high_index_material, low_index_material,
                    high_index_material, low_index_material,
                    high_index_material, low_index_material,
                    high_index_material
                ]

                thicknesses = [
                    high_thickness, low_thickness,
                    high_thickness, low_thickness,
                    high_thickness, low_thickness,
                    high_thickness
                ]

                # 不再反转顺序，因为我们已经修复了计算顺序问题
                # 膜层顺序：从入射介质到基材

                # 逐层添加
                for i, (material, thickness) in enumerate(
                        zip(materials, thicknesses)):
                    print(f"添加第{i+1}层: {material}, 厚度={thickness:.1f}nm")
                    self.add_layer()

                    if not self.layer_frames:
                        print(f"错误：无法添加第{i+1}层")
                        continue

                    layer = self.layer_frames[0]
                    layer.material_var.set(material)
                    layer.d_var.set(f"{thickness:.1f}")
                    layer.on_material_selected(None)

                # 基材设为玻璃
                self.substrate_var.set("Glass")
                self.substrate_thickness_var.set("1")
                self.substrate_unit_var.set("mm")
                if hasattr(self, 'on_substrate_selected'):
                    self.on_substrate_selected(None)

                # 出射介质设为空气
                self.exit_medium_var.set("Air")

            # 2. 宽带增透膜 (Broadband AR coating)
            elif "宽带增透膜" in preset_name or "broadband_ar" in preset_name_lower:
                print("正在应用宽带增透膜预设...")

                # 设置入射介质
                self.incident_medium_var.set("Air")

                # 三层宽带增透膜设计 - 使用固定的厚度值
                # 膜层顺序：从入射介质到基材
                materials = ["MgF2", "Al2O3", "HfO2"]
                thicknesses = [95.0, 45.0, 22.0]  # 优化的厚度，单位nm

                # 不再反转顺序，因为我们已经修复了计算顺序问题

                # 逐层添加
                for i, (material, thickness) in enumerate(
                        zip(materials, thicknesses)):
                    print(f"添加第{i+1}层: {material}, 厚度={thickness:.1f}nm")
                    self.add_layer()

                    if not self.layer_frames:
                        print(f"错误：无法添加第{i+1}层")
                        continue

                    layer = self.layer_frames[0]
                    layer.material_var.set(material)
                    layer.d_var.set(f"{thickness:.1f}")
                    layer.on_material_selected(None)

                # 基材设为玻璃
                self.substrate_var.set("Glass")
                self.substrate_thickness_var.set("1")
                self.substrate_unit_var.set("mm")
                if hasattr(self, 'on_substrate_selected'):
                    self.on_substrate_selected(None)

                # 出射介质设为空气
                self.exit_medium_var.set("Air")

            # 3. 窄带滤光片 (Narrow bandpass filter)
            elif "窄带滤光片" in preset_name or "narrow_band" in preset_name_lower:
                print("正在应用窄带滤光片预设...")

                # 设置入射介质
                self.incident_medium_var.set("Air")

                # 使用固定厚度值
                high_thickness = 57.0  # nm
                low_thickness = 94.0  # nm
                cavity_thickness = 188.0  # nm (半波长腔)

                # 高低折射率材料
                high_index_material = "TiO2"
                low_index_material = "SiO2"

                # 创建窄带滤光片 (HL)^3 2L (LH)^3 - 从入射介质到基材方向
                materials = [
                    # 第一个反射镜 (HL)^3
                    high_index_material, low_index_material,
                    high_index_material, low_index_material,
                    high_index_material, low_index_material,
                    # 半波长腔 2L
                    low_index_material,
                    # 第二个反射镜 (LH)^3
                    low_index_material, high_index_material,
                    low_index_material, high_index_material,
                    low_index_material, high_index_material
                ]

                thicknesses = [
                    # 第一个反射镜
                    high_thickness, low_thickness,
                    high_thickness, low_thickness,
                    high_thickness, low_thickness,
                    # 半波长腔
                    cavity_thickness,
                    # 第二个反射镜
                    low_thickness, high_thickness,
                    low_thickness, high_thickness,
                    low_thickness, high_thickness
                ]

                # 不再反转顺序，因为我们已经修复了计算顺序问题
                # 膜层顺序：从入射介质到基材

                # 逐层添加
                for i, (material, thickness) in enumerate(
                        zip(materials, thicknesses)):
                    print(f"添加第{i+1}层: {material}, 厚度={thickness:.1f}nm")
                    self.add_layer()

                    if not self.layer_frames:
                        print(f"错误：无法添加第{i+1}层")
                        continue

                    layer = self.layer_frames[0]
                    layer.material_var.set(material)
                    layer.d_var.set(f"{thickness:.1f}")
                    layer.on_material_selected(None)

                # 基材设为玻璃
                self.substrate_var.set("Glass")
                self.substrate_thickness_var.set("1")
                self.substrate_unit_var.set("mm")
                if hasattr(self, 'on_substrate_selected'):
                    self.on_substrate_selected(None)

                # 出射介质设为空气
                self.exit_medium_var.set("Air")

            # 4. 单层增透膜 (Single layer AR coating)
            elif "单层增透膜" in preset_name or "single_layer_ar" in preset_name_lower:
                print("正在应用单层增透膜预设...")

                # 设置入射介质
                self.incident_medium_var.set("Air")

                # 使用MgF2作为增透材料 - 从入射介质到基材方向
                material = "MgF2"

                # 使用固定厚度值
                thickness = 100.0  # nm

                # 添加单层
                self.add_layer()
                if self.layer_frames:
                    layer = self.layer_frames[0]
                    layer.material_var.set(material)
                    layer.d_var.set(f"{thickness:.1f}")
                    layer.on_material_selected(None)

                # 基材设为玻璃
                self.substrate_var.set("Glass")
                self.substrate_thickness_var.set("1")
                self.substrate_unit_var.set("mm")
                if hasattr(self, 'on_substrate_selected'):
                    self.on_substrate_selected(None)

                # 出射介质设为空气
                self.exit_medium_var.set("Air")

            # 5. 金属反射膜 (Metal reflector)
            elif "金属反射膜" in preset_name or "metal_reflector" in preset_name_lower:
                print("正在应用金属反射膜预设...")

                # 设置入射介质
                self.incident_medium_var.set("Air")

                # 使用铝作为反射材料 - 从入射介质到基材方向
                material = "Al"
                thickness = 100.0  # 100nm厚的铝膜

                # 添加金属层
                self.add_layer()
                if self.layer_frames:
                    layer = self.layer_frames[0]
                    layer.material_var.set(material)
                    layer.d_var.set(f"{thickness:.1f}")
                    layer.on_material_selected(None)

                # 基材设为玻璃
                self.substrate_var.set("Glass")
                self.substrate_thickness_var.set("1")
                self.substrate_unit_var.set("mm")
                if hasattr(self, 'on_substrate_selected'):
                    self.on_substrate_selected(None)

                # 出射介质设为空气
                self.exit_medium_var.set("Air")

            # 6. 低通高反膜 (Low-pass filter)
            elif "低通高反膜" in preset_name or "low_pass_filter" in preset_name_lower:
                print("正在应用低通高反膜预设...")

                # 设置入射介质
                self.incident_medium_var.set("Air")

                # 低通滤光片设计 - 高折射率材料在外侧，低折射率材料在内侧
                # 这样设计可以让低波长透过，高波长反射
                high_index_material = "TiO2"  # 高折射率材料
                low_index_material = "SiO2"   # 低折射率材料

                # 使用固定厚度值 - 针对低通特性优化
                high_thickness = 68.0  # nm
                low_thickness = 110.0  # nm

                # 创建低通滤光片结构 - 从入射介质到基材方向
                materials = []
                thicknesses = []

                # 添加8层交替膜系 (HLHLHLHL)
                for i in range(4):
                    materials.extend([high_index_material, low_index_material])
                    thicknesses.extend([high_thickness, low_thickness])

                # 逐层添加
                for i, (material, thickness) in enumerate(
                        zip(materials, thicknesses)):
                    print(f"添加第{i+1}层: {material}, 厚度={thickness:.1f}nm")
                    self.add_layer()

                    if not self.layer_frames:
                        print(f"错误：无法添加第{i+1}层")
                        continue

                    layer = self.layer_frames[0]
                    layer.material_var.set(material)
                    layer.d_var.set(f"{thickness:.1f}")
                    layer.on_material_selected(None)

            # 基材设为玻璃
            self.substrate_var.set("Glass")
            self.substrate_thickness_var.set("1")
            self.substrate_unit_var.set("mm")
            if hasattr(self, 'on_substrate_selected'):
                self.on_substrate_selected(None)

            # 出射介质设为空气
            self.exit_medium_var.set("Air")

            # 确保界面更新
            self.update_idletasks()
            self.update()

            print(
                f"已成功应用预设 '{preset_name}'，共添加{len(self.layer_frames)}层上表面膜层，{len(self.exit_layer_frames)}层下表面膜层")

        except Exception as e:
            print(f"应用预设时发生错误: {str(e)}")
            import traceback
            traceback.print_exc()
            messagebox.showerror("应用预设错误", f"应用预设时发生错误: {str(e)}")

    def manage_materials(self):
        """调用 MaterialManager 的 manage_materials 方法并刷新材料"""
        self.material_manager.manage_materials()
        self.wait_window(self.material_manager.window)
        
        # 重新加载材料数据
        load_materials_from_file("custom_materials.json")
        
        # 更新所有下拉框的选项
        material_names = get_material_names()
        try:
            self.incident_medium_combo['values'] = material_names
            self.substrate_combo['values'] = material_names
            self.exit_medium_combo['values'] = material_names
            for layer_frame in self.layer_frames + self.exit_layer_frames:
                layer_frame.material_combo['values'] = material_names
        except tk.TclError as e:
            print(f"更新下拉框时出错: {e}")
            self.recreate_combos()

    def update_combobox_values(self):
        material_names = get_material_names()
        try:
            # 合并所有下拉框引用
            combos = [self.incident_medium_combo, 
                    self.substrate_combo,
                    *[frame.material_combo for frame in self.layer_frames]]
            for combo in combos:
                combo['values'] = material_names
        except tk.TclError:
            self.recreate_combos()


if __name__ == "__main__":
    app = MultilayerFilmApp()
    app.mainloop()
