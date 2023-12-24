import matplotlib.pyplot as plt
from math import sqrt, log, log2, ceil

# Estimated time
# n = 1000000
# p = 2
# a = 0.001521 # network latency
# b = 1738499544 # network bandwidth
# clock_rate = 2.4 * 10**9 / 12 # 2.4 GHz
# time = sqrt(n) * log(log(sqrt(n))) / clock_rate + ((n * log(log(sqrt(n)))) / clock_rate) / 2*p + a * (ceil(log2(p)) + 1)
# print(f'time: {time}')

# calcualte karp flatt metric
processor = [1, 2, 4, 8, 16, 32, 64, 128]
speedup = [1,2,3.99,7.94,15.72,31,57.2,107.19]

karpFlatt = []
for i in range(1, len(processor)):
    karpFlatt.append((1/speedup[i] - 1/processor[i]) / (1 - 1/processor[i]))
print(karpFlatt)


# The performance of diagram
RealTime = [7.73257,3.86961,1.9391,0.97445,0.49181,0.2494,0.13519,0.07214]
estimateTime = [7.7325,3.8657,1.9384,0.97398,0.49284,0.2487,0.13481,0.07168]

plt.figure(figsize=(6, 4)) 
plt.plot(processor, RealTime, label='Real Execution Time', color='green')
plt.plot(processor, estimateTime, label='Estimate Execution Time', color='blue')
plt.xlabel('Processors')
plt.ylabel('Run time (second)')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.18))
plt.show()