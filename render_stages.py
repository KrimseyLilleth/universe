"""Render three key frames showing each stage of the involuted torus."""
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import math

W, H = 480, 480
PHI = 1.618
R_MAJOR = PHI
R_MINOR = 1.0

def involuted_torus_point(u, v, fold):
    theta = u * 2 * math.pi
    twist = fold * 0.5
    phi_angle = v * 2 * math.pi + theta * twist
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)
    envelope = 1.0 - fold * 0.45 * cos_t
    r_eff = R_MINOR * envelope
    center_deform = fold * R_MINOR * 0.5 * max(0, -cos_t)
    R_eff = R_MAJOR - center_deform
    cx = R_eff * cos_t
    cz = R_eff * sin_t
    x = cx + r_eff * math.cos(phi_angle) * cos_t
    y = r_eff * math.sin(phi_angle)
    z = cz + r_eff * math.cos(phi_angle) * sin_t
    return x, y, z

def project(x, y, z, rot_y, rot_x=0.35):
    cos_r = math.cos(rot_y)
    sin_r = math.sin(rot_y)
    x2 = x * cos_r - z * sin_r
    z2 = x * sin_r + z * cos_r
    cos_rx = math.cos(rot_x)
    sin_rx = math.sin(rot_x)
    y2 = y * cos_rx - z2 * sin_rx
    z3 = y * sin_rx + z2 * cos_rx
    d = 8.0
    scale = d / (d + z3)
    sx = W / 2 + x2 * scale * 100
    sy = H / 2 - y2 * scale * 100
    return sx, sy, z3, scale

def render_frame(fold, rot_y, show_particles=False, particle_time=0):
    img = Image.new('RGB', (W, H), (5, 5, 16))
    draw = ImageDraw.Draw(img)
    rng = np.random.RandomState(42)
    for _ in range(200):
        sx = rng.randint(0, W)
        sy = rng.randint(0, H)
        b = rng.randint(30, 80)
        draw.point((sx, sy), fill=(b, b, b + 20))

    nu, nv = 72, 36
    points = []
    for i in range(nu):
        for j in range(nv):
            u = i / nu
            v = j / nv
            x, y, z = involuted_torus_point(u, v, fold)
            sx, sy, depth, scale = project(x, y, z, rot_y)
            u_next = (i + 1) / nu
            v_next = (j + 1) / nv
            x2, y2, z2 = involuted_torus_point(u_next, v, fold)
            x3, y3, z3 = involuted_torus_point(u, v_next, fold)
            du = np.array([x2 - x, y2 - y, z2 - z])
            dv = np.array([x3 - x, y3 - y, z3 - z])
            normal = np.cross(du, dv)
            nl = np.linalg.norm(normal) + 1e-10
            normal = normal / nl
            vd = np.array([0, 0, -1])
            facing = abs(np.dot(normal, vd))
            fresnel = (1.0 - facing) ** 2
            t = u + v * 0.3
            base_r = int(78 + (255 - 78) * (math.sin(t * math.pi) * 0.5 + 0.5))
            base_g = int(205 + (209 - 205) * (math.sin(t * math.pi) * 0.5 + 0.5))
            base_b = int(196 + (102 - 196) * (math.sin(t * math.pi) * 0.5 + 0.5))
            base_r = int(base_r + (255 - base_r) * fold * fresnel * 0.5)
            base_g = int(base_g * (1 - fold * fresnel * 0.3))
            base_b = int(base_b * (1 - fold * fresnel * 0.3))
            light = 0.3 + 0.5 * facing + 0.4 * fresnel
            r = min(255, int(base_r * light))
            g = min(255, int(base_g * light))
            b = min(255, int(base_b * light))
            points.append((sx, sy, depth, r, g, b, scale))

    points.sort(key=lambda p: -p[2])
    for sx, sy, depth, r, g, b, scale in points:
        if 0 <= sx < W and 0 <= sy < H:
            size = max(1, int(2.5 * scale))
            draw.ellipse([sx - size, sy - size, sx + size, sy + size], fill=(r, g, b))

    if show_particles:
        rng2 = np.random.RandomState(123)
        for i in range(600):
            pu = (rng2.random() + particle_time * 0.08) % 2.0
            pv = (rng2.random() + particle_time * 0.02) % 1.0
            x, y, z = involuted_torus_point(pu * 0.5, pv, fold)
            sx, sy, depth, scale = project(x, y, z, rot_y)
            if 0 <= sx < W and 0 <= sy < H:
                brightness = 0.5 + 0.5 * math.sin(pu * math.pi + particle_time)
                pr = min(255, int(200 * brightness + 55))
                pg = min(255, int(255 * brightness))
                pb = min(255, int(180 * brightness + 75))
                ps = max(1, int(2 * scale))
                draw.ellipse([sx - ps - 2, sy - ps - 2, sx + ps + 2, sy + ps + 2],
                           fill=(pr // 4, pg // 4, pb // 4))
                draw.ellipse([sx - ps, sy - ps, sx + ps, sy + ps], fill=(pr, pg, pb))
    return img

def add_labels(img, title, subtitle, bottom_text=""):
    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf", 22)
        font_sm = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf", 14)
        font_bot = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf", 13)
    except:
        font = ImageFont.load_default()
        font_sm = font
        font_bot = font
    draw.rectangle([12, 12, 460, 95], fill=(5, 5, 16))
    draw.text((18, 18), title, fill=(255, 209, 102), font=font)
    draw.text((18, 48), subtitle, fill=(170, 170, 170), font=font_sm)
    draw.text((18, 70), f"R/r = \u03c6 = {PHI}", fill=(255, 209, 102), font=font_sm)
    if bottom_text:
        draw.rectangle([12, H - 45, W - 12, H - 12], fill=(5, 5, 16))
        draw.text((18, H - 40), bottom_text, fill=(78, 205, 196), font=font_bot)
    return img

# Stage 1
print("Rendering Stage 1...")
img1 = render_frame(fold=0, rot_y=0.6)
img1 = add_labels(img1, "Stage 1: Normal Torus",
    "A surface with clear inside and outside",
    "The starting shape — a standard torus with R/r = \u03c6")
img1.save("/home/user/universe/stage1_torus.png")

# Stage 2
print("Rendering Stage 2...")
img2 = render_frame(fold=0.55, rot_y=1.0)
img2 = add_labels(img2, "Stage 2: Folding Through",
    "The surface pushes through its own wall",
    "One side shrinks, the other expands through the center")
img2.save("/home/user/universe/stage2_fold.png")

# Stage 3
print("Rendering Stage 3...")
img3 = render_frame(fold=1.0, rot_y=1.4, show_particles=True, particle_time=3.0)
img3 = add_labels(img3, "Stage 3: Involuted Flow",
    "No inside or outside — continuous surface",
    "Particles need TWO full loops to return to start")
img3.save("/home/user/universe/stage3_flow.png")

print("Done!")
