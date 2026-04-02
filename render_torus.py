"""Render an animated GIF of the involuted torus in three stages."""
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import math

W, H = 480, 480
FRAMES_STAGE1 = 20   # normal torus rotating
FRAMES_STAGE2 = 35   # folding animation
FRAMES_STAGE3 = 25   # involuted with particles
TOTAL = FRAMES_STAGE1 + FRAMES_STAGE2 + FRAMES_STAGE3

PHI = 1.618
R_MAJOR = PHI
R_MINOR = 1.0

def involuted_torus_point(u, v, fold):
    """Compute a point on the involuted torus surface."""
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
    """Project 3D point to 2D with rotation and perspective."""
    # Rotate around Y
    cos_r = math.cos(rot_y)
    sin_r = math.sin(rot_y)
    x2 = x * cos_r - z * sin_r
    z2 = x * sin_r + z * cos_r

    # Rotate around X (tilt)
    cos_rx = math.cos(rot_x)
    sin_rx = math.sin(rot_x)
    y2 = y * cos_rx - z2 * sin_rx
    z3 = y * sin_rx + z2 * cos_rx

    # Perspective
    d = 8.0
    scale = d / (d + z3)
    sx = W / 2 + x2 * scale * 100
    sy = H / 2 - y2 * scale * 100

    return sx, sy, z3, scale

def render_frame(fold, rot_y, show_particles=False, particle_time=0):
    """Render one frame of the torus."""
    img = Image.new('RGB', (W, H), (5, 5, 16))
    draw = ImageDraw.Draw(img)

    # Draw background stars
    rng = np.random.RandomState(42)
    for _ in range(200):
        sx = rng.randint(0, W)
        sy = rng.randint(0, H)
        b = rng.randint(30, 80)
        draw.point((sx, sy), fill=(b, b, b + 20))

    # Generate mesh points
    nu, nv = 64, 32
    points = []

    for i in range(nu):
        for j in range(nv):
            u = i / nu
            v = j / nv
            x, y, z = involuted_torus_point(u, v, fold)
            sx, sy, depth, scale = project(x, y, z, rot_y)

            # Color based on position and fold
            # Fresnel-like: edges brighter
            u_next = (i + 1) / nu
            v_next = (j + 1) / nv
            x2, y2, z2 = involuted_torus_point(u_next, v, fold)
            x3, y3, z3 = involuted_torus_point(u, v_next, fold)

            # Normal approximation
            du = np.array([x2 - x, y2 - y, z2 - z])
            dv = np.array([x3 - x, y3 - y, z3 - z])
            normal = np.cross(du, dv)
            nl = np.linalg.norm(normal) + 1e-10
            normal = normal / nl

            # View direction (towards camera)
            vd = np.array([0, 0, -1])
            facing = abs(np.dot(normal, vd))
            fresnel = 1.0 - facing
            fresnel = fresnel ** 2

            # Base color: cyan-gold gradient
            t = u + v * 0.3
            base_r = int(78 + (255 - 78) * (math.sin(t * math.pi) * 0.5 + 0.5))
            base_g = int(205 + (209 - 205) * (math.sin(t * math.pi) * 0.5 + 0.5))
            base_b = int(196 + (102 - 196) * (math.sin(t * math.pi) * 0.5 + 0.5))

            # Mix with red during fold
            base_r = int(base_r + (255 - base_r) * fold * fresnel * 0.5)
            base_g = int(base_g * (1 - fold * fresnel * 0.3))
            base_b = int(base_b * (1 - fold * fresnel * 0.3))

            # Lighting
            light = 0.3 + 0.5 * facing + 0.4 * fresnel
            r = min(255, int(base_r * light))
            g = min(255, int(base_g * light))
            b = min(255, int(base_b * light))

            points.append((sx, sy, depth, r, g, b, scale))

    # Sort by depth (painter's algorithm)
    points.sort(key=lambda p: -p[2])

    # Draw points
    for sx, sy, depth, r, g, b, scale in points:
        if 0 <= sx < W and 0 <= sy < H:
            size = max(1, int(2.5 * scale))
            draw.ellipse([sx - size, sy - size, sx + size, sy + size], fill=(r, g, b))

    # Draw flow particles in stage 3
    if show_particles:
        rng2 = np.random.RandomState(123)
        for i in range(600):
            pu = (rng2.random() + particle_time * 0.08) % 2.0
            pv = (rng2.random() + particle_time * 0.02) % 1.0
            x, y, z = involuted_torus_point(pu * 0.5, pv, fold)

            # Offset slightly from surface
            sx, sy, depth, scale = project(x, y, z, rot_y)

            if 0 <= sx < W and 0 <= sy < H:
                # Bright glowing particle
                brightness = 0.5 + 0.5 * math.sin(pu * math.pi + particle_time)
                pr = min(255, int(200 * brightness + 55))
                pg = min(255, int(255 * brightness))
                pb = min(255, int(180 * brightness + 75))
                ps = max(1, int(2 * scale))
                # Glow
                draw.ellipse([sx - ps - 2, sy - ps - 2, sx + ps + 2, sy + ps + 2],
                           fill=(pr // 4, pg // 4, pb // 4))
                draw.ellipse([sx - ps, sy - ps, sx + ps, sy + ps], fill=(pr, pg, pb))

    return img

def add_labels(img, text, subtitle=""):
    """Add text overlay to frame."""
    draw = ImageDraw.Draw(img)

    # Title
    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf", 20)
        font_sm = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf", 14)
    except:
        font = ImageFont.load_default()
        font_sm = font

    # Dark background for text
    draw.rectangle([15, 15, 400, 90], fill=(5, 5, 16, 200))

    # Gold title
    draw.text((20, 20), text, fill=(255, 209, 102), font=font)

    if subtitle:
        draw.text((20, 48), subtitle, fill=(170, 170, 170), font=font_sm)

    # Phi label
    draw.text((20, 68), f"R/r = φ = {PHI}", fill=(255, 209, 102), font=font_sm)

    return img

print("Rendering frames...")
frames = []

# Stage 1: Normal torus rotating
for i in range(FRAMES_STAGE1):
    rot = i / FRAMES_STAGE1 * math.pi * 0.5 + 0.3
    img = render_frame(fold=0, rot_y=rot)
    img = add_labels(img, "Stage 1: Normal Torus", "A surface with inside and outside")
    frames.append(img)
    if i % 10 == 0:
        print(f"  Stage 1: {i}/{FRAMES_STAGE1}")

# Stage 2: Folding
for i in range(FRAMES_STAGE2):
    t = i / FRAMES_STAGE2
    fold = t  # linear 0->1
    rot = 0.3 + math.pi * 0.5 + t * math.pi * 0.3
    img = render_frame(fold=fold, rot_y=rot)
    pct = int(fold * 100)
    img = add_labels(img, "Stage 2: Involution", f"Folding through itself... {pct}%")
    frames.append(img)
    if i % 10 == 0:
        print(f"  Stage 2: {i}/{FRAMES_STAGE2}")

# Stage 3: Involuted with particles
for i in range(FRAMES_STAGE3):
    t = i / FRAMES_STAGE3
    rot = 0.3 + math.pi * 0.8 + t * math.pi * 0.5
    img = render_frame(fold=1.0, rot_y=rot, show_particles=True, particle_time=t * 10)
    img = add_labels(img, "Stage 3: Involuted Flow", "Two loops to return — no inside or outside")
    frames.append(img)
    if i % 10 == 0:
        print(f"  Stage 3: {i}/{FRAMES_STAGE3}")

print("Saving GIF...")
output_path = "/home/user/universe/involuted-torus.gif"
frames[0].save(
    output_path,
    save_all=True,
    append_images=frames[1:],
    duration=100,  # ms per frame
    loop=0
)

print(f"Done! Saved to {output_path}")
print(f"Total frames: {len(frames)}")
import os
size_mb = os.path.getsize(output_path) / (1024 * 1024)
print(f"File size: {size_mb:.1f} MB")
