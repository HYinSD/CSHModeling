# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 17:15:38 2023

@author: 41569
"""

import tkinter as tk
import numpy as np
import random
# Modification item
Cao = 2

#
X_time = 2
Y_time = 2
Z_time = 2
root = tk.Tk()
coordinate = np.loadtxt('inputcoord.txt')


# Setting window titles, creating labels and text boxes
root.title("Get user input of Ca/Si ratio")
root.geometry("450x100")
label = tk.Label(root, text="Please enter a calcium-silicon ratio greater than 1.2 and less than 2.3.：")
label.pack()
entry = tk.Entry(root)
entry.pack()
# Interlayer calcium charge
if Cao ==2:
    mask = (coordinate[:, 3] == 1.36)
    coordinate[mask, 3] += 0.64
    mask = (coordinate[:, 3] == 1.05)
    coordinate[mask, 3] += 0.31
    
####################
if Cao ==1.36:
    mask = (coordinate[:, 3] == 1.36)
    coordinate[mask, 3] += 0.64
    mask = (coordinate[:, 3] == 1.05)
    coordinate[mask, 3] += 0.31    
# Cell Expansion
Expansion = np.array([[0,0,0,0],
                  [0,0,0,0]])
for i1 in range(1,X_time):
    for i2 in range(1, Y_time):
        for i3 in range(1, Z_time):
            x = 22.317317677222220*(i1-1)        
            y = 22.177195157222220*(i2-1)
            z = 22.770053081222223*(i3-1)
            add = np.array([x,y,z,0])
            both = coordinate + add
            Expansion = np.append(Expansion,both,axis = 0)
Expansion = Expansion[~np.all(Expansion == 0, axis=1)]
# Calculate the number of target aggregates
def submit():
    value = entry.get()
    try:
        CS = float(value)
        def calculate_values(CS):
            if 1.2 <= CS <= 1.25:
                Q0, Q1, Q2 = 0.0179211469534050, 0.387096774193548, 0.594982078853047
            elif 1.25 < CS <= 1.35:
                Q0, Q1, Q2 = 0.0107526881720429, 0.573476702508961, 0.415770609318996
            elif 1.35 < CS <= 1.45:
                Q0, Q1, Q2 = 0.0143369175627240, 0.691756272401434, 0.293906810035842
            elif 1.45 < CS <= 1.55:
                Q0, Q1, Q2 = 0.0107526881720430, 0.759856630824373, 0.229390681003584
            elif 1.55 < CS <= 1.65:
                Q0, Q1, Q2 = 0.0107526881720430, 0.767025089605735, 0.222222222222222
            elif 1.65 < CS <= 1.75:
                Q0, Q1, Q2 = 0.00716845878136198, 0.788530465949821, 0.204301075268817
            elif 1.75 < CS <= 1.85:
                Q0, Q1, Q2 = 0.0143369175627240, 0.784946236559140, 0.200716845878136
            elif 1.85 < CS <= 1.95:
                Q0, Q1, Q2 = 0.0143369175627241, 0.817204301075269, 0.168458781362007
            elif 1.95 < CS <= 2.05:
                Q0, Q1, Q2 = 0.0107526881720430, 0.845878136200717, 0.143369175627240
            elif 2.05 < CS <= 2.15:
                Q0, Q1, Q2 = 0.0143369175627240, 0.870967741935484, 0.114695340501792
            elif 2.15 < CS <= 2.25:
                Q0, Q1, Q2 = 0.0107526881720430, 0.888888888888889, 0.100358422939068
            elif 2.25 < CS <= 2.3:
                Q0, Q1, Q2 = 0.0107526881720429, 0.910394265232975, 0.0788530465949821
            else:
                return None
            N1 = round((144 / CS) * Q0)
            N5 = round(((144 / CS) * Q2) / 3)
            N2 = round((((144 / CS) * Q1) - 2 * N5) / 2)

            return N1, N2, N5
        N1, N2, N5 = calculate_values(CS)  
        N1, N2, N5 = calculate_values(CS) 
        root.destroy()  # Close window
        # Atomic variable naming 
        cah = Expansion == 1.36  
        ob = Expansion == -1.050
        obos = Expansion == -1.1808
        cao = Expansion == Cao
        st = Expansion == 2.100

        # Determine the number of molecules 
        numcah = np.sum(Expansion == 1.36)
        numst = np.sum(Expansion == 2.100)
        numobos = np.sum(Expansion == Expansion.min())
        numcao = np.sum(cao)
        numob = np.size(Expansion)/4-numcah-numst-numobos-numcao
        #########################################################################
        datc = np.array([[0,0,0,0],
                         [0,0,0,0]])
        yos = np.where(Expansion[:,3]==2.1)
        dstc = Expansion[yos]
        data = dstc[np.argsort(dstc[:,1])]
        for i4 in range(4*(X_time-1)):
            for i6 in range(Z_time-1):
                aond1 = data[:,0] >= 1.39481019810000-2.28195+5.5793*i4
                aond2 = data[:,0] <=2.41025017410000+2.28195+5.5793*i4
                aond3 = data[:,2] >=2.26563054010000-2.28195+22.770053081222223*i6
                aond4 = data[:,2] <=4.20107969010000+2.28195+22.770053081222223*i6
                aond = (aond1 & aond2 & aond3 & aond4)
                dstc1 = data[aond]

                bond1 = data[:,0] >= 1.39481019810000-2.28195+5.5793*i4
                bond2 = data[:,0] <=2.32098072410000+2.28195+5.5793*i4
                bond3 = data[:,2] >=7.41164828010000-2.28195+22.770053081222223*i6
                bond4 = data[:,2] <=9.48371737010000+2.28195+22.770053081222223*i6
                bond = (bond1 & bond2 & bond3 & bond4)
                dstc2 = data[bond]

                cond1 = data[:,0] >=3.16904048410000-2.28195+5.5793*i4
                cond2 = data[:,0] <=4.18448045910000+2.28195+5.5793*i4
                cond3 = data[:,2] >=13.6506555400000-2.28195+22.770053081222223*i6
                cond4 = data[:,2] <=15.5861046900000+2.28195+22.770053081222223*i6
                cond = (cond1 & cond2 & cond3 & cond4)
                dstc3 = data[cond]
                dond1 = data[:,0] >=3.25830993310000-2.28195+5.5793*i4
                dond2 = data[:,0] <=4.18448045810000+2.28195+5.5793*i4
                dond3 = data[:,2] >=18.7966732800000-2.28195+22.770053081222223*i6
                dond4 = data[:,2] <=20.8459723800000+2.28195+22.770053081222223*i6
                dond = (dond1 & dond2 & dond3 & dond4)
                dstc4 = data[dond]
                dstc = np.vstack((dstc1,dstc2,dstc3,dstc4))
                datc = np.append(datc,dstc,axis = 0)
        datc = datc[~np.all(datc == 0, axis=1)]
        # Number of clusters based on distance
        distance = 3.15
        count = 0
        data = datc.reshape((16, 9,-1))
        while   True:
            q1 = 0
            q2 = 0
            q5 = 0
            count = 0
            q1matrices=np.array([[0,0,0,0],
                           [0,0,0,0]])
            q2matrices=np.array([[0,0,0,0],
                           [0,0,0,0]])
            q5matrices=np.array([[0,0,0,0],
                           [0,0,0,0]])
            for i in range(16):
                n1=np.array([[2,3,4],
                             [1,4,7],
                             [2,4,7],
                             [3,6,9],
                             [1,4,7],
                             [6,7,8],
                             [6,8,9], 
                             [2,3,4],
                             [1,2,4],
                             [2,5,8],
                             [3,5,7],
                             [3,6,8],
                             [2,4,7],])
                n1 = n1 - np.array([1,1,1])
                n2=np.array([[1,4],
                             [6,8], 
                             [2,4],
                             [6,7],
                             [6,9],
                             [3,4],
                             [1,4],
                             [2,8],
                             [1,7],
                             [3,9]])
                n2 = n2 - np.array([1,1])
                chosen1 = random.sample(range(len(n1)), 1)
                chosen1 = n1[chosen1]
                chosen2 = random.sample(range(len(n2)), 1)
                chosen2 = n2[chosen2]
                chosen = random.choice([chosen1, chosen2])
                remove = np.delete(data[i],chosen, axis=0)
                count = 0
                for n in range(remove.shape[0]-1):
                    gap = np.linalg.norm(remove[n+1,:] - remove[n,:])
                    if gap < distance:
                        count += 1
                    else:
                         if count >= 0: 
                             if count+1 == 1:
                                 q1 =q1+1
                                 q1matrices = np.vstack((q1matrices,remove[n,:]))
                             elif count+1 == 2:
                                 q2 =q2+1   
                                 q2matrices = np.vstack((q2matrices,remove[n-1:n+1,:]))
                             elif count+1 ==5:
                                 q5 =q5+1     
                                 q5matrices = np.vstack((q5matrices,remove[n-4:n+1,:]))
                         count = 0
                if count >= 0:
                    if count+1 == 1:
                        q1 =q1+1
                        q1matrices = np.vstack((q1matrices,remove[n+1,:]))
                    elif count+1 == 2:
                        q2 =q2+1   
                        q2matrices = np.vstack((q2matrices,remove[n:n+2,:]))
                    elif count+1 ==5:
                        q5 =q5+1
                        q5matrices = np.vstack((q5matrices,remove[n-3:n+2,:]))
                if q5 == N5 and q2== N2 and q1 >= N1: 
                    break 
            if q5 == N5 and q2 == N2 and q1 >= N1:
                break
        q1matrices = q1matrices[~np.all(q1matrices == 0, axis=1)]
        idx = np.random.choice(q1matrices.shape[0], N1, replace=False)
        q1matrices = q1matrices[idx,:]
        q2matrices = q2matrices[~np.all(q2matrices == 0, axis=1)] 
        q5matrices = q5matrices[~np.all(q5matrices == 0, axis=1)]
        st = np.vstack((q1matrices,q2matrices,q5matrices))

        #The oxygen atom of the polymer determines
        o1 = np.where(Expansion[:,3]==-1.05)
        o2 =np.where(Expansion[:,3]==-1.1808)
        ob= Expansion[o1]
        obos = Expansion[o2]
        ob = ob[np.argsort(ob[:,1])]
        obos = obos[np.argsort(obos[:,1])]
        O = np.append(ob,obos,axis = 0)
        number = 0

        xyz = O[:,0:3]
        xy = st[:,0:3]
        xyz1 = q1matrices[:,0:3]
        xyz2 =q2matrices[:,0:3]
        xyz5 = q5matrices[:,0:3]

        o1matrices = np.array([[0,0,0,0],
                       [0,0,0,0]])
        o2matrices =np.array([[0,0,0,0],
                       [0,0,0,0]])
        o5matrices =np.array([[0,0,0,0],
                       [0,0,0,0]])
        for x in range(len(q1matrices)):
            for y in range(len(O)):
                distance = np.linalg.norm(xyz[y,:] - xyz1[x,:])
                if distance < 2:
                    o1matrices = np.vstack((o1matrices,O[y,:]))
        o1matrices = o1matrices[~np.all(o1matrices == 0, axis=1)]
        o1matrices[:,3] = -1.1808
        for z in range(len(q2matrices)):
             for l in range(len(O)):
                distance = np.linalg.norm(xyz[l,:] - xyz2[z,:])
                if distance < 2:
                    o2matrices = np.vstack((o2matrices,O[l,:]))
        o2matrices = o2matrices[~np.all(o2matrices == 0, axis=1)]
        o2matrices = np.unique(o2matrices,axis=0)
        for z in range(len(q5matrices)):
            for l in range(len(O)):
                distance = np.linalg.norm(xyz[l,:] - xyz5[z,:])
                if distance < 2:
                    o5matrices = np.vstack((o5matrices,O[l,:]))
        o5matrices = o5matrices[~np.all(o5matrices == 0, axis=1)]
        o5matrices = np.unique(o5matrices,axis=0)
        o5matrices = np.vstack((o2matrices,o5matrices))
        o5j = o5matrices[:,0:3]
        for i in range(len(o5j)):
            for u in range(len(xy)):
                distance = np.linalg.norm(o5j[i,:]-xy[u,:])
                if distance < 2:
                    number = number +1
                    if number == 2:
                        o5matrices[i,3] = -1.05
                    else:
                        o5matrices[i,3] = -1.1808
            number = 0
        omatrices = np.vstack((o1matrices,o5matrices))
        #Statistical charge
        #########################################################################
        #ca
        cs = np.where(Expansion[:,3]==1.36)
        cs1 =np.where(Expansion[:,3]== Cao)
        ca = Expansion[cs]
        ca1 = Expansion[cs1]
        ca = np.append(ca1,ca,axis = 0)
        num_rows_to_delete = 40
        rows_to_delete = np.random.choice(ca.shape[0], num_rows_to_delete, replace=False)
        ca = np.delete(ca, rows_to_delete, axis=0)
        cad = np.sum(ca[:,3],axis = 0)
        #st

        std = np.sum(st[:,3],axis = 0)
        #o
        od = np.sum(omatrices[:,3],axis = 0)
        #electric charge difference
        need = cad + std +od
        #########################################################################
        # equilibrium charge
        # Remaining free oxygen
        mask = np.isin(O, omatrices).all(axis=1)
        oho = O[~mask]
        oho[:,3]=-0.95
        # equilibrium charge
        oh = oho.copy()

        def generate_random_point_nearby(original_point):
            random_vector = np.random.randn(3)
            norm = np.linalg.norm(random_vector)
            if norm != 0:
                random_vector /= norm
            else:
                random_vector = np.random.randn(3)
            new_point = original_point[:3] + random_vector  
            final_distance = np.linalg.norm(new_point - original_point[:3])  
            if not np.isclose(final_distance, 1.0):
                new_point = original_point[:3] + (new_point - original_point[:3]) / final_distance  
            new_point_with_425 = np.append(new_point, 0.425)  
            return new_point_with_425


        new_oh = []
        for i in range(len(oho)):
            point = oho[i]  # Use oho here instead of O
            new_point = generate_random_point_nearby(point)
            new_point_reshaped = new_point.reshape(1, -1)
    
            new_oh.append(np.append(point[:3], -0.95)) 
            new_oh.append(new_point_reshaped[0])  # new point based on oho

        oh = np.array(new_oh)


        n = need/0.525
        n = round(2 * n)
        supplementalO = oh[:n]




        # All atoms
        qo = np.vstack((st,omatrices,ca))
        num_rows, num_cols = qo.shape
        new_col = np.arange(num_rows)
        qo = np.insert(qo, 0, new_col+1, axis=1)
        add = np.ones((num_rows,1))
        qo = np.insert(qo, 1, add[1], axis=1)
        qo = qo[:,[0,1,5,2,3,4]]
        col = qo[:, 2]
        new_col = np.where(col == 1.36, 1,
                  np.where(col == 2.1, 2,
                  np.where(col == -1.05, 3,
                  np.where(col == -1.1808, 4,
                  np.where(col == -0.95, 6,
                  np.where(col == 0.425, 7,
                  np.where(col == Cao, 5,0)))))))
        qo = np.insert(qo, 2, new_col, axis=1)


        end = np.vstack((st,omatrices,ca,supplementalO))




        num_rows, num_cols = end.shape
        new_col = np.arange(num_rows)
        end = np.insert(end, 0, new_col+1, axis=1)
        add = np.ones((num_rows,1))
        end = np.insert(end, 1, add[1], axis=1)
        end = end[:,[0,1,5,2,3,4]]
        col = end[:, 2]


        #########################################################################
        new_col = np.where(col == 1.36, 1,
                  np.where(col == 2.1, 2,
                  np.where(col == -1.05, 3,
                  np.where(col == -1.1808, 4,
                  np.where(col == -0.95, 6,
                  np.where(col == 0.425, 7,
                  np.where(col == Cao, 5,0)))))))
        end = np.insert(end, 2, new_col, axis=1)
        #########################################################################
        count = np.sum(end[:, 2] == 6)

        col3_eq_6 = end[end[:, 2] == 6, 0]  # 第三列=6时，取第一列数据
        col3_eq_7 = end[end[:, 2] == 7, 0]  # 第三列=7时，取第一列数据

# 确保两个数组长度相同（如果不同，补齐较短的）
        max_len = max(len(col3_eq_6), len(col3_eq_7))
        col3_eq_6 = np.pad(col3_eq_6, (0, max_len - len(col3_eq_6)), constant_values=0)  # 用 0 补齐
        col3_eq_7 = np.pad(col3_eq_7, (0, max_len - len(col3_eq_7)), constant_values=0)

# 生成最终矩阵
        matrix = np.column_stack((np.arange(1, max_len + 1), np.ones(max_len), col3_eq_6, col3_eq_7))
        print(q1matrices)
        print(q2matrices)
        print(q5matrices)
        # File Export
        np.savetxt('output.data', end,fmt='%7.0f %6.0f %3.0f %9.6f %15.9f %15.9f %15.9f',delimiter=' ')
        with open('output.data', 'r') as f:
            content = f.read()
        with open('output.data', 'w') as f:
            f.write("LAMMPS data file, CSHModeling / CGCMM for Hamid\n")
            
            
            f.write(str(len(end)) + ' atoms\n')
            f.write(f"{count} bonds \n")
            f.write( """ 
            0 angles

           7 atom types
           1 bond types
           0 angle types

             0.000009703    22.317309703 xlo xhi
             0.022721230    22.199921230 ylo yhi
             0.116710684    22.886810684 zlo zhi

        Masses

           1  40.080000 # cah
           2  28.085500 # st
           3  15.999400 # ob
           4  15.999400 # obos
           5  40.080000 # cao
           6  15.999400
           7  1
           
        Atoms # full
        \n""")
            f.write(''.join(content) )
            f.close()
        with open("output.data", "a") as f:
            f.write("\n") 
            f.write("Bonds\n\n")  # 写入 "Bonds"
            np.savetxt(f, matrix, fmt="%d")  # 以整数格式写入矩阵
    except ValueError:
        print("Error")
# Create a submit button
button = tk.Button(root, text="OK", command=submit)
button.pack(pady=20)

root.mainloop()
