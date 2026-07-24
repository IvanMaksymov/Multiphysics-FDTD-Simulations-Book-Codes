import serial
import matplotlib.pyplot as plt
import numpy as np
import time

PORT = "/dev/ttyACM1"
BAUD = 115200
KE = 200

ser = serial.Serial(PORT, BAUD, timeout=1)

time.sleep(2)
ser.reset_input_buffer()

x = np.arange(KE)

plt.ion()

fig, ax = plt.subplots()
line, = ax.plot(x, np.zeros(KE))

ax.set_xlim(0, KE)
ax.set_ylim(-0.2, 1.1)

plt.show()

while True:

    raw = ser.readline().decode("ascii", errors="ignore").strip()

    if not raw:
        continue

    try:
        y = np.array([float(v) for v in raw.split(",")])

        if len(y) == KE:
            line.set_ydata(y)
            fig.canvas.draw_idle()
            fig.canvas.flush_events()

    except ValueError:
        continue
