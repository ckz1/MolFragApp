import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import colorsys

# Define the rgbcolor class
class rgbcolor:
    def __init__(self, initlist):
        excluded = [[0.12, 0.22]]  # exclude yellow hues from the colorwheel
        excluded.sort(key=lambda x: x[0])
        temp1 = [[min(1., max(0., el[0])), max(0., min(1., el[1]))] for el in excluded]
        temp2 = [[0., 0.]]
        for el in temp1:
            if el[0] >= temp2[-1][1]:
                temp2.append(el)
            else:
                temp2[-1][1] = el[1]
        self.excluded = temp2[1:]
        self.initlist = [max(0, el) for el in initlist]
        self.n = sum(el > 0 for el in self.initlist)
        self.m = len(initlist)
        self.a = 1. - sum(el[1] - el[0] for el in self.excluded)
        self.startlist = [0.] * self.m
        self.incrlist = [0.] * self.m
        for i in range(1, self.m):
            self.startlist[i] = self.startlist[i - 1] + (self.a / self.n if self.initlist[i - 1] > 0 else 0)
        for i in range(self.m):
            if self.initlist[i] > 0:
                self.incrlist[i] = self.a / self.n / self.initlist[i]

    def rgb_to_hex(self, rgb):
        return '#{:02x}{:02x}{:02x}'.format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))

    def hexcolor(self, index, el):
        if not (1 <= index <= self.m and 1 <= el <= self.initlist[index - 1]):
            return '#FFFFFF'
        hue = self.startlist[index - 1] + self.incrlist[index - 1] * (el - 1)
        for start, end in self.excluded:
            if hue > start:
                hue += end - start
        return self.rgb_to_hex(colorsys.hsv_to_rgb(hue, 1, 1))

# 用户可配置参数
plot_range = (0, 16)
bin_count = 30
kde_bandwidth = 0.2
opacity = 0.3

# 加载数据
df = pd.read_csv('2024-11-15T08-01_export.csv')

color_generator = rgbcolor([1] * len(df['state'].unique()))
state_colors = [color_generator.hexcolor(i + 1, 1) for i in range(len(df['state'].unique()))]

# 创建图形并调整子图比例，设置sharex=True以共享x轴
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 16), gridspec_kw={'height_ratios': [0.7, 0.3]}, sharex=True)

y_kde_sum = np.zeros(1000)
x_range = np.linspace(plot_range[0], plot_range[1], 1000)

for state, color in zip(df['state'].unique(), state_colors):
    subset = df[df['state'] == state]['delta_energy']
    hist, edges = np.histogram(subset, bins=np.linspace(plot_range[0], plot_range[1], bin_count))
    ax2.bar(edges[:-1], hist, width=np.diff(edges)[0], alpha=opacity, color=color, label=state)
    
    kde = gaussian_kde(subset, bw_method=kde_bandwidth)
    y_kde = kde.evaluate(x_range)
    ax1.fill_between(x_range, 0, y_kde, alpha=opacity, color=color, label=state)
    y_kde_sum += y_kde

# 归一化并绘制总和曲线
y_kde_sum_normalized = y_kde_sum / y_kde_sum.max()
ax1.plot(x_range, y_kde_sum_normalized, color='black', label='Normalized Total Sum', linewidth=2)

# 设置图表的整体属性
ax1.set_ylabel('Intensity (arb. units)')
ax1.set_ylim(0, 1.1)
ax1.legend()

ax2.set_ylabel('Count')
ax2.set_ylim(0, 11.5)
ax2.set_xlabel('KER (eV)')
fig.subplots_adjust(hspace=0)  # 移除子图之间的空白

plt.show()