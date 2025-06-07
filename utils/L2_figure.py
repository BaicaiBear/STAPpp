import numpy as np
import matplotlib.pyplot as plt

# 问题 (4): 计算L2误差（线性单元）
h_list = [1, 0.5, 0.25, 0.125]
errors_linear = [2.62437e-01, 1.28634e-01, 5.95555e-02, 3.26347e-02]

# 绘制双对数图
log_h = np.log(h_list)
log_err_linear = np.log(errors_linear)
plt.figure()
plt.plot(log_h, log_err_linear, 'o-', label='Linear Element')
plt.xlabel('log(h)')
plt.ylabel('log(L2 Error)')
plt.title('S4R Element Error Convergence')
slope, _ = np.polyfit(log_h, log_err_linear, 1)
plt.plot(log_h, slope*log_h, '--', label=f'Slope={slope:.2f}')
plt.legend()
plt.grid(True)
plt.show()