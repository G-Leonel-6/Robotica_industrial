import numpy as np
import scipy as sc
import roboticstoolbox as rtb
import spatialmath as sm
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from matplotlib.animation import FuncAnimation
import time
import imageio
import imageio.v3 as iio
import imageio_ffmpeg as ffmpeg
import os
from PIL import Image

class irb140_clase(rtb.DHRobot):
    def __init__(self):
        super().__init__([
            rtb.RevoluteDH(alpha=-np.pi/2,a=0.07, qlim=[-np.pi,np.pi]),
            rtb.RevoluteDH(a=0.36,offset=-np.pi/2, qlim=[-np.pi,np.pi]),
            rtb.RevoluteDH(alpha=np.pi/2,offset=np.pi, qlim=[-np.pi,np.pi]),
            rtb.RevoluteDH(d=0.38, alpha=-np.pi/2, qlim=[-np.pi,np.pi]),
            rtb.RevoluteDH(alpha=np.pi/2, qlim=[-np.pi,np.pi]),
            rtb.RevoluteDH()
        ], name="IRB140")

    def get_config(self,q):
        g1 = np.sign(self.links[3].d * np.sin(q[1]+self.links[1].offset +q[2]+self.links[2].offset) + self.links[1].a * np.cos(q[1]+self.links[1].offset) + self.links[0].a)
        g2 = np.sign(np.cos(q[2]+self.links[2].offset))
        g3 = np.sign(np.sin(q[4]+self.links[4].offset))
        return np.array([g1,g2,g3])

    def get_joint_bounds(self):
        return np.array([link.qlim if link.qlim is not None else [-np.pi, np.pi] for link in self.links])

    def solve_q1(self, g1, q2):
        d3, o1, o2, a1, a0 = self.links[3].d, self.links[1].offset, self.links[2].offset, self.links[1].a, self.links[0].a
        def S(theta): return d3 * np.sin(theta + q2 + o2) + a1 * np.cos(theta) + a0
        sol = fsolve(S, 0.0)[0]
        theta = sol + np.pi/6 if g1 * S(sol + np.pi/6) > 0 else sol - np.pi/6
        return theta - o1

    def guess_seed_from_conf(self, conf,POSE ,q_ref=None):
        q0 = np.zeros(6)

        px, py, pz = POSE.t
        g1 = conf[0]
        if px >= 0 and py >= 0: 
            q0[0] = np.pi/4 if g1 == 1 else -3*np.pi/4
        elif px >= 0 and py < 0:  
            q0[0] = -np.pi/4 if g1 == 1 else 3*np.pi/4
        elif px < 0 and py >= 0:  
            q0[0] = 3*np.pi/4 if g1 == 1 else -np.pi/4
        else:  
            q0[0] = -3*np.pi/4 if g1 == 1 else np.pi/4

        q0[2] = np.pi/2 + conf[1] * np.pi/3
        q0[4] = conf[2] * np.pi/4
        q0[1] = self.solve_q1(conf[0], q0[2])
        return 0.5*q0 + 0.5*q_ref if q_ref is not None else q0

    def wrap_joints(self, q):
        q_wrapped = np.copy(q)
        bounds = self.get_joint_bounds()
        for i in range(len(q)):
            q_min, q_max = bounds[i]
            span = q_max - q_min
            q_wrapped[i] = (q[i] - q_min) % span + q_min
        return q_wrapped

    def damped_pinv(self, J, damping=1e-4):
        U, S, Vh = np.linalg.svd(J)
        S_inv = S / (S**2 + damping**2)
        return Vh.T @ np.diag(S_inv) @ U.T

    def config_match(self, q, conf):
        return np.all(self.get_config(q) == conf)

    def ikine_recta_con_trayectoria(self, POSE_deseada, conf, q_inicial=None,
                                    vd=0.01, Ts=0.001, tol=1e-3, max_iter=2000,
                                    q_actual=None):
        q_i = q_inicial if q_inicial is not None else (
            q_actual.copy() if q_actual is not None else self.guess_seed_from_conf(conf, POSE_deseada)
        )
        trayectoria = [q_i.copy()]
        errores = []

        for i in range(max_iter):
            POSE_i = self.fkine(q_i)
            e_vec = sm.base.tr2delta(POSE_i.A, POSE_deseada.A)
            error_norm = np.linalg.norm(e_vec)
            errores.append(error_norm)

            if error_norm < tol:
                return (q_i, 0, np.array(trayectoria), errores) if self.config_match(q_i, conf) else (q_i, -2, np.array(trayectoria), errores)

            delta_pose_i = (e_vec / error_norm) * vd * Ts
            J = self.jacobe(q_i)
            dq_i = self.damped_pinv(J) @ delta_pose_i
            q_i = self.wrap_joints(q_i + dq_i)
            trayectoria.append(q_i)

        return q_i, -1, np.array(trayectoria), errores

    def ikine_recta_con_trayectoria_multi(self, POSE_deseada, conf, q_actual=None,
                                          vd=0.01, Ts=0.001, tol=1e-3, max_iter=2000,
                                          n_reintentos=6):
        for intento in range(n_reintentos):
            if intento == 0 and q_actual is not None:
                q_seed = q_actual
            elif intento == 1:
                q_seed = self.guess_seed_from_conf(conf, POSE_deseada)
            else:
                perturb = (np.random.rand(6)-0.5)*0.2
                q_seed = self.guess_seed_from_conf(conf, POSE_deseada) + perturb

            q, status, trayectoria, errores = self.ikine_recta_con_trayectoria(
                POSE_deseada, conf, q_inicial=q_seed, vd=vd, Ts=Ts, tol=tol, max_iter=max_iter
            )
            if status == 0:
                return q, status, trayectoria, errores
        return q, status, trayectoria, errores
    

robot = irb140_clase()
n_acierto = 0
n_iter = 1

for i in range(n_iter):
    q_deseado = (np.random.rand(6)-0.5)*2*np.pi
    POSE = robot.fkine(q_deseado)
    conf = robot.get_config(q_deseado)
    q, success, history, _ = robot.ikine_recta_con_trayectoria_multi(
        POSE, conf=conf, vd=1, Ts=0.001, max_iter=8000, n_reintentos=10
    )
    if success == 0:
        n_acierto += 1
    else:
        print(f"❌ ERROR:\nDeseado: {q_deseado}\nConfig: {conf}\nAlcanzado: {q}")

print("✅ Cantidad de corridas:", n_iter)
print("🎯 Cantidad de aciertos:", n_acierto)

def getframe_safe(fig):
    # Forzar render completo
    fig.canvas.draw_idle()
    plt.pause(0.1)  # pequeña pausa para permitir actualización del backend

    # Ahora obtener datos del canvas
    buf = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    w, h = fig.canvas.get_width_height()
    expected_size = w * h * 3

    # Si no coincide, ajustamos el tamaño automáticamente
    if buf.size != expected_size:
        # intentamos estimar altura
        h = int(buf.size / (w * 3))
    buf = buf.reshape(h, w, 3)
    buf = np.flipud(buf)
    return Image.fromarray(buf)

frames = []
for i, q in enumerate(history):
    robot.plot(q, block=False)
 # aseguramos render completo
    frame = getframe_safe(plt.gcf())
    plt.savefig(f"frames/frame_{i:04d}.png")
    plt.clf()


"""
env = rtb.backends.PyPlot.PyPlot()
env.launch()

env.add(robot)
frames = []



output_dir = "frames_animacion"
os.makedirs(output_dir, exist_ok=True)

for i, frame in enumerate(frames):
    filename = os.path.join(output_dir, f"frame_{i:04d}.png")
    frame.save(filename)
    print(f"Guardado: {filename}")

print(f"\n✅ Se guardaron {len(frames)} frames en la carpeta: {output_dir}")




def getframe_fixed(fig):
    fig.canvas.draw()
    # Obtener bytes RGB
    buf = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    # Calcular tamaño real a partir de los ejes del canvas
    w, h = fig.canvas.get_width_height()
    expected_size = w * h * 3

    # Si no coincide, ajustamos el tamaño automáticamente
    if buf.size != expected_size:
        # intentamos estimar altura
        h = int(buf.size / (w * 3))
        print(f"[getframe_fixed] Ajustando tamaño detectado: ({w}, {h})")

    # Dar forma y voltear
    buf = buf.reshape(h, w, 3)
    buf = np.flipud(buf)
    return Image.fromarray(buf)

for q in history:
    robot.q = q
    env.step(0.01)
    frame = getframe_fixed(plt.gcf())
    frames.append(frame)

output_dir = "frames_animacion"
os.makedirs(output_dir, exist_ok=True)

for i, frame in enumerate(frames):
    filename = os.path.join(output_dir, f"frame_{i:04d}.png")
    frame.save(filename)
    print(f"Guardado: {filename}")

print(f"\n✅ Se guardaron {len(frames)} frames en la carpeta: {output_dir}")

"""