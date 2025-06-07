from pylab import *
import numpy as np
from scipy.interpolate import griddata
from subprocess import Popen, PIPE

def split_vector(v,m,n):
    ans = []
    r = []
    for j in range(n):
        for i in range(m-1):
            pos = i*n + j
            r.append(v[pos])
        ans.append(r)
        r = []
    return ans

def read_data():
    lines = None
    with open("output.txt","r") as file:
        lines = file.readlines()
    r_i,r_e,m,n,iso,ninst = lines[0].split()
    r_i = float(r_i)
    r_e = float(r_e)
    m = int(m)
    n = int(n)
    iso = float(iso)
    ninst = int(ninst)
    len_entry = (m-1)*n + n - 1
    data = []
    temp_vector = []
    isos = []
    right = 0
    for i in range(ninst):
        line = lines
        left = right + 1
        right = left + len_entry
        for j in range(left,right + 1 - n):
            line = lines[j].rstrip("\n")
            temp_vector.append(float(line))
        for j in range(right + 1 - n,right + 1):
            line = lines[j].rstrip("\n")
            isos.append(int(line))
        data.append((temp_vector,isos))
        temp_vector = []
        isos = []
    new_data = []
    for p in data:
        v = p[0]
        iso_arr = p[1]
        v = split_vector(v,m,n)
        iso_arr = [r_i + x*(r_e - r_i)/m for x in iso_arr]
        new_data.append((v,iso_arr))
    print(new_data)
    for i,p in enumerate(new_data):
        t,r_iso = p
        graph_temperature(n,m-1,r_i,r_e,t,experiment=i)
        graph_iso(r_iso,n,r_i,r_e,iso,ticks=5,experiment=i)

def graph_iso(r_isos:list,theta_partitions,r_i,r_e,temperature=500,ticks = 5,experiment = 0):
     # Solo para que la isoterma se "pegue" bien al dar la vuelta
    isos = r_isos.copy()
    isos.append(isos[0])
    theta = np.linspace(0, 2*np.pi, theta_partitions + 1)
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(theta, isos)
    
    ax.set_rticks([r_i + (r_e - r_i) * x / ticks for x in range(ticks + 1)])
    ax.grid(True)
    ax.set_title(f"Isoterma {temperature} (Experimento {experiment})")
    plt.show()

def graph_temperature(theta_partitions,r_partitions,r_i,r_e,temperatures,experiment=0):
    temperatures.append(temperatures[0])
    theta = np.tile(np.linspace(0, 2*np.pi, theta_partitions + 1), (r_partitions, 1)).transpose()
    r = np.tile(np.linspace(r_i, r_e, r_partitions), (theta_partitions + 1, 1))

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ticks = [100*i for i in range(1,16)]
    ticknames = [str(tick) for tick in ticks]
    grafico = ax.pcolor(theta, r, temperatures, cmap='jet')
    ax.set_title(f'Alto horno (Experimento {experiment})')

    ax.set_yticklabels(ticknames)

    fig.colorbar(grafico,ticks=ticks)
    plt.show()
def graph_data(r_partitions,theta_partitions,points,t_values):
    #create 5000 Random points distributed within the circle radius 100
    


    #Some function to generate values for these points, 
    #this could be values = np.random.rand(number_points)
    r_i = points[0][0]
    r_e = points[len(points) - 1][0]
    points = np.array(points)
    values = np.array(t_values)
    print(values)

    #now we create a grid of values, interpolated from our random sample above
    theta = np.linspace(0.0, 2*np.pi, theta_partitions)
    r = np.linspace(r_i,r_e, r_partitions)
    print(r)
    grid_r, grid_theta = np.meshgrid(r, theta)
    data = griddata(points, values, (grid_r, grid_theta), method='linear')

    #Create a polar projection
    ax1 = plt.subplot(projection="polar")
    ax1.pcolormesh(theta, r, data.T,cmap="hot")
    plt.plot()
    plt.show()

def main():
    r_partitions,theta_partitions,points,t_values = read_data()
    graph_data(r_partitions,theta_partitions,points,t_values)



if __name__ == "__main__":
    pass