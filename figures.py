import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------------------------ #
# Figure 1
# ------------------------------------------------------------------------------------------------ #

df = pd.read_csv('q3_data.csv')

dataFrames = []
dataFrames.append(df.loc[df['runtime'] == 1000])
dataFrames.append(df.loc[df['runtime'] == 2000])
dataFrames.append(df.loc[df['runtime'] == 3000])

for i in range(3):
    plt.plot(dataFrames[i]["utilization"], dataFrames[i]["average_size"])
    
plt.legend(["Runtime = 1000s", "Runtime = 2000s", "Runtime = 3000s"])
plt.title("Average Queue Size vs Average Utilization")
plt.xlabel("Utilization")
plt.ylabel("Average Queue Size (packets)")

plt.savefig("figure1.jpg")
plt.clf()

# ------------------------------------------------------------------------------------------------ #
# Figure 2
# ------------------------------------------------------------------------------------------------ #

for i in range(3):
    plt.plot(dataFrames[i]["utilization"], dataFrames[i]["idle_proportion"])
    
plt.legend(["Runtime = 1000s", "Runtime = 2000s", "Runtime = 3000s"])
plt.title("Idle Proportion of Queue vs Average Utilization")
plt.xlabel("Utilization")
plt.ylabel("Idle Proportion")

plt.savefig("figure2.jpg")
plt.clf()

# ------------------------------------------------------------------------------------------------ #
# Figure 3
# ------------------------------------------------------------------------------------------------ #

dataFrames = []
temp = df.loc[df['runtime'] == 3000]
dataFrames.append(temp.loc[temp['capacity'] == 10])
dataFrames.append(temp.loc[temp['capacity'] == 25])
dataFrames.append(temp.loc[temp['capacity'] == 50])
dataFrames.append(temp.loc[temp['capacity'] == 75])

for i in range(4):
    plt.plot(dataFrames[i]["utilization"], dataFrames[i]["average_size"])

plt.legend([
    "Queue Capacity = 10 packets", 
    "Queue Capacity = 25 packets", 
    "Queue Capacity = 50 packets",
    "Queue Capacity = 75 packets"])

plt.title("Average Queue Size vs Queue Utilization")
plt.xlabel("Utilization")
plt.ylabel("Average Queue Size (packets)")

plt.savefig("figure3.jpg")
plt.clf()

# ------------------------------------------------------------------------------------------------ #
# Figure 4
# ------------------------------------------------------------------------------------------------ #

for i in range(4):
    plt.plot(dataFrames[i]["utilization"], dataFrames[i]["drop_proportion"])

plt.legend([
    "Queue Capacity = 10 packets", 
    "Queue Capacity = 25 packets", 
    "Queue Capacity = 50 packets",
    "Queue Capacity = 75 packets"])

plt.title("Proportion of Received Packet Dropped vs Queue Utilization")
plt.xlabel("Utilization")
plt.ylabel("Drop Proportion")

plt.savefig("figure4.jpg")
plt.clf()

# ------------------------------------------------------------------------------------------------ #
# Figure 5 and 6
# ------------------------------------------------------------------------------------------------ #

df = pd.read_csv('q6_data.csv')

dataFrames = []
temp = df.loc[df['capacity'] == 25]
dataFrames.append(temp.loc[temp['runtime'] == 1000])
dataFrames.append(temp.loc[temp['runtime'] == 2000])
dataFrames.append(temp.loc[temp['runtime'] == 3000])

for i in range(3):
    plt.plot(dataFrames[i]["utilization"], dataFrames[i]["average_size"])

plt.legend(["Runtime = 1000s", "Runtime = 2000s", "Runtime = 3000s"])
plt.title("Average Queue Size vs Average Utilization")
plt.xlabel("Utilization")
plt.ylabel("Average Queue Size (packets)")

plt.savefig("figure5.jpg")
plt.clf()

for i in range(3):
    plt.plot(dataFrames[i]["utilization"], dataFrames[i]["drop_proportion"])
    
plt.legend(["Runtime = 1000s", "Runtime = 2000s", "Runtime = 3000s"])
plt.title("Queue Drop Proportion vs Average Utilization")
plt.xlabel("Utilization")
plt.ylabel("Drop Proportion")

plt.savefig("figure6.jpg")
plt.clf()
