# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 23:27:17 2025

@author: 41569
"""

import tkinter as tk
from tkinter import messagebox
import subprocess

def run_script(script_name):

    try:
        subprocess.run(["python", script_name])
    except Exception as e:
        messagebox.showerror("error", f"Error occurred while executing the script: {e}")

def on_select(option):

    if option == "1":
        run_script("Obos-Based Charge Compensation.py")
    elif option == "2":
        run_script("Hydroxyl-Based Charge Compensation.py")  
    else:
        messagebox.showerror("error", "Invalid selection")

def create_gui():

    window = tk.Tk()

    label = tk.Label(window, text="Please select the method for charge compensation.:")
    label.pack(pady=10)

    button1 = tk.Button(window, text="Obos-Based Charge Compensation", command=lambda: on_select("1"))
    button1.pack(pady=5)
    
    button2 = tk.Button(window, text="Hydroxyl-Based Charge Compensation", command=lambda: on_select("2"))
    button2.pack(pady=5)

    window.mainloop()

if __name__ == "__main__":
    create_gui()
