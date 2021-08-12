import csv
import matplotlib.pyplot as plt

ts = []
xs = []
with open("First_particle_20Hz.csv") as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        t = float(row[0])
        x = float(row[1])
        ts.append(t)
        xs.append(x)

plt.figure(figsize=(6.4,4.8))
plt.plot(ts, xs, color='blue')
plt.xlabel("Time (s)")
plt.ylabel("Displacement (mm)")
plt.title("Position of the First Particle - Experiment", fontsize='small')
plt.grid(color='skyblue', linestyle=':', linewidth=0.5)
plt.tight_layout()
plt.show()

ts = []
xs = []
with open("Second_particle_20Hz.csv") as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        t = float(row[0])
        x = float(row[1])
        ts.append(t)
        xs.append(x)

plt.figure(figsize=(6.4,4.8))
plt.plot(ts, xs, color='blue')
plt.xlabel("Time (s)")
plt.ylabel("Displacement (mm)")
plt.title("Position of the Second Particle - Experiment", fontsize='small')
plt.grid(color='skyblue', linestyle=':', linewidth=0.5)
plt.tight_layout()
plt.show()

ts = []
xs = []
with open("Third_particle_20Hz.csv") as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        t = float(row[0])
        x = float(row[1])
        ts.append(t)
        xs.append(x)

plt.figure(figsize=(6.4,4.8))
plt.plot(ts, xs, color='blue')
plt.xlabel("Time (s)")
plt.ylabel("Displacement (mm)")
plt.title("Position of the Third Particle - Experiment", fontsize='small')
plt.grid(color='skyblue', linestyle=':', linewidth=0.5)
plt.tight_layout()
plt.show()