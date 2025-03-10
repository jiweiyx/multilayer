import tkinter as tk
from tkinter import ttk, messagebox
import json
from material_database import MATERIALS, load_materials_from_file, get_material_names, get_material_n_k_data, get_material_description
import os
import sys

def resource_path(relative_path):
    """获取数据文件的绝对路径"""
    try:
        # PyInstaller创建的临时文件夹中运行
        base_path = sys._MEIPASS
    except Exception:
        # Use the directory where manager.py is located
        base_path = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base_path, relative_path)

class MaterialManager:
    def __init__(self, parent=None):
        self.parent = parent
        self.window = None
    def on_close(self):
        self.parent.update_combobox_values()  # 调用主窗口更新方法
        self.window.destroy()
        self.window = None

    def manage_materials(self):
        if not self.window or not self.window.winfo_exists():  # 双重校验
            self.window = tk.Toplevel(self.parent)
            self.window.protocol("WM_DELETE_WINDOW", self.on_close)
        self.window.title("Manage Materials")
        self.window.geometry("700x450")

        tk.Label(self.window, text="Select Material:").grid(row=0, column=0, padx=5, pady=5)
        self.material_var = tk.StringVar()
        material_dropdown = ttk.Combobox(self.window, textvariable=self.material_var, values=get_material_names())
        material_dropdown.grid(row=0, column=1, padx=5, pady=5)
        material_dropdown.current(0)

        load_button = tk.Button(self.window, text="Load", command=lambda: self.load_material(self.material_var.get()))
        load_button.grid(row=0, column=2, padx=5, pady=5)

        new_button = tk.Button(self.window, text="New", command=self.new_material)
        new_button.grid(row=0, column=3, padx=5, pady=5)

        tk.Label(self.window, text="Description:").grid(row=1, column=0, padx=5, pady=5)
        self.desc_var = tk.StringVar()
        tk.Entry(self.window, textvariable=self.desc_var, width=40).grid(row=1, column=1, columnspan=3, padx=5, pady=5)

        tk.Label(self.window, text="Data Format:").grid(row=1, column=4, padx=5, pady=5)
        self.format_var = tk.StringVar(value="table")
        tk.Entry(self.window, textvariable=self.format_var, width=10, state='readonly').grid(row=1, column=5, padx=5, pady=5)

        self.tree = ttk.Treeview(self.window, columns=("Wavelength", "n", "k"), show="headings", height=12)
        self.tree.heading("Wavelength", text="Wavelength (nm)")
        self.tree.heading("n", text="Refractive Index (n)")
        self.tree.heading("k", text="Extinction Coefficient (k)")
        self.tree.column("Wavelength", width=150)
        self.tree.column("n", width=150)
        self.tree.column("k", width=150)
        self.tree.grid(row=2, column=0, columnspan=6, padx=5, pady=5)

        self.tree.bind("<Double-1>", self.edit_cell)

        add_row_button = tk.Button(self.window, text="Add Row", command=self.add_row)
        add_row_button.grid(row=3, column=0, padx=5, pady=5)

        delete_row_button = tk.Button(self.window, text="Delete Row", command=self.delete_row)
        delete_row_button.grid(row=3, column=1, padx=5, pady=5)

        save_button = tk.Button(self.window, text="Save Data", command=self.save_data)
        save_button.grid(row=3, column=2, columnspan=2, padx=5, pady=5)

        delete_button = tk.Button(self.window, text="Delete Data", command=self.delete_data)
        delete_button.grid(row=3, column=4, columnspan=2, padx=5, pady=5)

        self.load_material(material_dropdown.get())

        if self.parent:
            self.window.transient(self.parent)
            self.window.grab_set()
        else:
            self.window.mainloop()

    def load_material(self, material_name):
        """加载材料数据到 UI"""
        if material_name in MATERIALS:
            self.current_material = material_name
            material_data = get_material_n_k_data(material_name)
            self.desc_var.set(get_material_description(material_name))
            self.format_var.set(MATERIALS[material_name]['data'])

            for item in self.tree.get_children():
                self.tree.delete(item)

            if material_data and 'wavelength' in material_data:
                for i in range(len(material_data['wavelength'])):
                    n_value = material_data['n'][i] if material_data['n'][i] is not None else ""
                    k_value = material_data['k'][i] if material_data['k'][i] is not None else ""
                    self.tree.insert("", "end", values=(
                        material_data['wavelength'][i],
                        n_value,
                        k_value
                    ))
            else:
                print(f"材料 {material_name} 的数据为空或格式错误")
        else:
            print(f"材料 {material_name} 不存在")

    def new_material(self):
        """创建新材料条目，并自动生成默认波长行"""
        new_name = f"NewMaterial_{len(MATERIALS)}"
        MATERIALS[new_name] = {
            'description': '',
            'data': 'table',
            'values': {}
        }
        print(f"新建空白材料: {new_name}")

        material_dropdown = self.window.winfo_children()[1]
        material_dropdown['values'] = get_material_names()
        material_dropdown.set(new_name)
        
        self.current_material = new_name
        self.desc_var.set('')
        for item in self.tree.get_children():
            self.tree.delete(item)

        default_wavelengths = [400, 450, 500, 550, 600, 650, 700, 750, 800]
        for wl in default_wavelengths:
            self.tree.insert("", "end", values=(wl, "", ""))

    def add_row(self):
        """添加新行到表格"""
        children = self.tree.get_children()
        if children:
            last_wl = float(self.tree.item(children[-1], 'values')[0])
            new_wl = last_wl + 50
        else:
            new_wl = 400
        self.tree.insert("", "end", values=(new_wl, "", ""))

    def delete_row(self):
        """删除选中的表格行"""
        selected_item = self.tree.selection()
        if not selected_item:
            messagebox.showwarning("Warning", "Please select a row to delete!")
            return
        self.tree.delete(selected_item[0])

    def edit_cell(self, event):
        """双击编辑表格中的任意列，确保输入框在原地"""
        item = self.tree.identify_row(event.y)
        if not item:
            return

        column = self.tree.identify_column(event.x)
        col_index = int(column[1:]) - 1

        values = self.tree.item(item, 'values')
        x, y, width, height = self.tree.bbox(item, column)
        if not (x and y and width and height):
            return

        entry = tk.Entry(self.tree, width=15)
        entry.place(x=x, y=y, width=width, height=height)
        entry.insert(0, values[col_index] if values[col_index] else "")
        entry.focus_set()

        def save_edit(event=None):
            try:
                new_value = float(entry.get()) if entry.get() else ""
                updated_values = list(values)
                updated_values[col_index] = new_value
                self.tree.item(item, values=tuple(updated_values))
                entry.destroy()
            except ValueError:
                messagebox.showerror("Error", "Please enter a valid number!")
                entry.focus_set()

        entry.bind("<Return>", save_edit)
        entry.bind("<FocusOut>", save_edit)

    def save_data(self):
        """保存当前材料数据，支持名称修改"""
        # 安全获取下拉框引用 [^2]
        material_dropdown = next(child for child in self.window.winfo_children() 
                            if isinstance(child, ttk.Combobox))
        new_name = material_dropdown.get().strip()
        
        # 基础校验
        if not hasattr(self, 'current_material'):
            messagebox.showerror("Error", "No material selected!")
            return
        if not new_name:
            messagebox.showerror("Error", "Material name cannot be empty!")
            return
        if new_name == 'Air':
            messagebox.showerror("Error", "Cannot use reserved name 'Air'!")
            return
        
        original_name = self.current_material
        
        # 处理重命名
        if new_name != original_name:
            if new_name in MATERIALS:  # [^1]
                messagebox.showerror("Error", f"Material '{new_name}' already exists!")
                return
            # 执行重命名操作
            MATERIALS[new_name] = MATERIALS.pop(original_name)
            self.current_material = new_name  # 更新当前材料名称

        # 数据校验
        values = {'wavelength': [], 'n': [], 'k': []}
        has_error = False
        for item in self.tree.get_children():
            vals = self.tree.item(item, 'values')
            wl, n, k = vals

            try:
                wl_float = float(wl) if wl else None
                n_float = float(n) if n else None
                k_float = float(k) if k else 0.0

                if not all([wl_float, n_float]):
                    messagebox.showwarning("Warning", f"Skipping invalid row: {wl}")
                    has_error = True
                    continue

                values['wavelength'].append(wl_float)
                values['n'].append(n_float)
                values['k'].append(k_float)
            except ValueError as e:
                messagebox.showerror("Error", f"Invalid data: {wl}, {e}")
                has_error = True

        if has_error or not values['wavelength']:
            messagebox.showerror("Error", "保存失败：存在无效数据")
            return

        # 更新数据并保存
        MATERIALS[self.current_material].update({ 
            'description': self.desc_var.get(),
            'values': values
        })
        
        # 持久化到文件
        with open(resource_path("custom_materials.json"), 'w', encoding='utf-8') as f:
            json.dump(MATERIALS, f, ensure_ascii=False, indent=4) 

        # 更新界面
        updated_names = get_material_names()
        material_dropdown['values'] = updated_names
        material_dropdown.set(new_name)
        
        # 同步主界面 [^4]
        if hasattr(self, 'parent'):
            for combo in [self.parent.incident_medium_combo, 
                        self.parent.substrate_combo]:
                combo['values'] = updated_names
                if combo.get() == original_name:
                    combo.set(new_name)

        messagebox.showinfo("Success", "材料数据保存成功！")


    def delete_data(self):
        """删除当前材料"""
        if not hasattr(self, 'current_material') or self.current_material == 'Air':
            messagebox.showerror("Error", "Cannot delete base material 'Air' or no material selected!")
            return

        if messagebox.askyesno("Confirm", f"Delete {self.current_material}?"):
            # 删除内存数据并保存到文件
            del MATERIALS[self.current_material]
            full_path = resource_path("custom_materials.json")
            with open(full_path, 'w', encoding='utf-8') as f:
                json.dump(MATERIALS, f, ensure_ascii=False, indent=4)
            
            # 关键：重新加载材料数据到内存 [^4]
            load_materials_from_file(full_path)  # 重新初始化MATERIALS
            
            messagebox.showinfo("Success", "Material data deleted!")

            # 安全获取当前窗口的下拉框
            material_dropdown = next(child for child in self.window.winfo_children() 
                                if isinstance(child, ttk.Combobox))
            
            # 更新所有相关下拉框（包括主界面）[^4]
            remaining_materials = get_material_names()
            # 当前窗口下拉框更新
            material_dropdown['values'] = remaining_materials
            # 主界面下拉框同步更新
            if hasattr(self, 'parent'):
                self.parent.incident_medium_combo['values'] = remaining_materials
                self.parent.substrate_combo['values'] = remaining_materials
            
            # 处理当前选择状态
            if remaining_materials:
                default_material = remaining_materials[0]
                material_dropdown.set(default_material)
                self.load_material(default_material)
            else:
                material_dropdown.set('')
                self.tree.delete(*self.tree.get_children())


if __name__ == "__main__":
    manager = MaterialManager()
    manager.manage_materials()