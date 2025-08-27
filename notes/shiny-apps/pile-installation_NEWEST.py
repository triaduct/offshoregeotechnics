import numpy as np
import matplotlib.pyplot as plt


"""
Downward direction is negative

At the moment, no elasticity is included
"""
#%%
class SoilLayer:
    def __init__(self, top, bottom, Qs_max, damping, shaft_quake):
        self.top = top
        self.bottom = bottom
        self.Qs_max = Qs_max                # ultimate shaft resistance [N]
        self.damping = damping      # Smith damping coefficient[s/m]
        self.shaft_quake = shaft_quake          # [m]

    def contains(self, z):
        return self.bottom <= z < self.top


class SoilProfile:
    def __init__(self, layers):
        self.layers = layers  # List of SoilLayer objects

    def get_layer_at_depth(self, z):
        for layer in self.layers:
            if layer.contains(z):
                return layer
        return SoilLayer(top=1e6, bottom=-1e6, Qs_max=0.0, damping=0.0, shaft_quake=1e-6)  # No soil present (e.g. above ground surface)

class SpringDashpot:
    def __init__(self, damping, max_force, quake):
        self.k = max_force / quake if quake > 0 else 1e6  # avoid div by zero
        self.Js = damping
        self.max_force = max_force
        self.quake = quake
        self.u_perm = 0.0
        self.R_total = 0.0

    def compute_force(self, u_rel, v_rel):
        delta_u = u_rel - self.u_perm
        sign = np.sign(delta_u)
        abs_du = abs(delta_u)

        if abs_du <= self.quake:
            R_static = self.k * delta_u
        else:
            R_static = self.max_force * sign
            self.u_perm += (abs_du - self.quake) * sign

        R_damp = self.Js * abs(R_static) * v_rel
        self.R_total = R_static + R_damp
        return self.R_total


class HAMMER:
    def __init__(self, mass, efficiency, E_rated, cushion_stiffness, cushion_damping, cushion_quake):
        self.mass = mass
        self.efficiency = efficiency
        self.E_rated = E_rated
        self.k_cushion = SpringDashpot(cushion_damping, E_rated / cushion_quake, cushion_quake)

    def impact_velocity(self, pile_mass):
        return -np.sqrt(2 * self.efficiency * self.E_rated / (pile_mass + self.mass))

    def reset(self, v0=0.0):
        self.k_cushion.u_perm = 0.0
        self.u = 0.0
        self.v = v0
    
class PILE:
    def __init__(self, length, diameter, bottom_elevation, n_increments, 
                 youngs_modulus, damping_ratio,
                 material_density,
                 soil_profile: SoilProfile):
        # Constants
        self.length = length
        self.diameter = diameter
        self.area = np.pi * (diameter / 2)**2
        self.n = n_increments
        self.dz = length / n_increments     # Length of each increment
        self.total_mass = self.area * length * material_density 
        self.mass_per_increment = self.total_mass / self.n
        self.m = np.ones(self.n) * self.mass_per_increment  # Mass vector 

        # Structural stiffness and damping. Assuming perfectly uniform pile.
        self.k_axial = (youngs_modulus * self.area) / self.dz   # N/m
        self.c_axial = 2 * damping_ratio * np.sqrt(self.k_axial * self.mass_per_increment)

        # # Shaft spring-dashpot elements (REMOVED: defined later)
        # self.k = [SpringDashpot(
        #             damping=0.16,           # s/m
        #             max_force=200e3,       # ultimate shaft resistance [N]
        #             quake=0.002            # 2 mm quake
        #         ) for _ in range(self.n)]
        
        # Base spring-dashpot element 
        self.base = SpringDashpot(
                    damping=0.16,           # s/m
                    max_force=1e6*0.5*20*np.pi*(0.5/2)**2,       # ultimate base resistance
                    quake=0.1*self.diameter             # 0.1D quake
                )

        # Dynamic state placeholders (Topmost is first)
        self.u = np.zeros(self.n)   # Displacement vector
        self.v = np.zeros(self.n)   # Velocity vector
        self.a = np.zeros(self.n)   # Acceleration vector

        # Blow tracking
        self.blow_depths = []

        # Elevation of each segment centerline (topmost is first)
        self.segment_elevations = [
            bottom_elevation + (i + 0.5) * self.dz for i in range(self.n)
        ]
        self.segment_elevations = self.segment_elevations[::-1] # Flip around
        self.assign_soil_springs(soil_profile)

        self.history = {
            "blow": [],
            "time": [],
            "top_disp": [],
            "top_vel": [],
            "top_acc": [],
            "tip_disp": [],
            "tip_vel": [],
        }
        
        self.blow_histories = []

    def assign_soil_springs(self, soil_profile: SoilProfile):
        self.k = []
        for i, z in enumerate(self.segment_elevations):
            layer = soil_profile.get_layer_at_depth(z)
            if layer is not None:
                # Spring active in this layer
                self.k.append(SpringDashpot(
                    damping=layer.damping,
                    max_force=layer.Qs_max * self.dz * self.diameter,  # perimeter * dz * fs
                    quake=layer.shaft_quake
                ))
            else:
                # No resistance (e.g. above ground)
                self.k.append(SpringDashpot(
                    damping=0.0,
                    max_force=0.0,
                    shaft_quake=1e-6  # small quake to avoid division by zero
                ))

    def compute_internal_forces(self):
        F_internal = np.zeros(self.n)
        for i in range(1, self.n):
            du = self.u[i] - self.u[i-1]
            dv = self.v[i] - self.v[i-1]
            F = (self.k_axial * du + self.c_axial * dv)      # Pile structural spring model
            F_internal[i] += F
            F_internal[i-1] -= F
        return F_internal

    def compute_shaft_resistances(self):
        return np.array([
            self.k[i].compute_force(self.u[i], self.v[i]) for i in range(self.n)
        ]) #[N]
 
    def compute_base_resistances(self):
        return self.base.compute_force(self.u[-1], self.v[-1]) # [N]

    def conservation_of_momentum(self, dt, hammer=None):
        F_int = self.compute_internal_forces()
        R_shaft = self.compute_shaft_resistances()
        R_base = self.compute_base_resistances()
        F_shifted = np.append(F_int[1:], R_base)
        self.a = (-F_int - self.m * 9.81 + R_shaft + F_shifted) / self.m
        self.v += self.a * dt
        self.u += self.v * dt

        if hammer:
            du = hammer.u - self.u[0]
            dv = hammer.v - self.v[0]
            Fcushion = hammer.k_cushion.compute_force(du, dv)
            a_hammer = (-Fcushion - hammer.mass * 9.81) / hammer.mass
            hammer.v += a_hammer * dt
            hammer.u += hammer.v * dt

            self.a[0] += Fcushion / self.m[0]
            self.v[0] += (Fcushion / self.m[0]) * dt
            self.u[0] += self.v[0] * dt

        t_now = self.history["time"][-1] + dt if self.history["time"] else 0
        self.history["time"].append(t_now)
        self.history["top_disp"].append(self.u[0])
        self.history["top_vel"].append(self.v[0])
        self.history["top_acc"].append(self.a[0])
        self.history["tip_disp"].append(self.u[-1])
        self.history["tip_vel"].append(self.v[-1])


    def simulate(self, hammer, n_blows, dt, t_per_blow, rest_time):
        for blow in range(n_blows):
            blow_time, blow_force = [], []

            for elem in self.k:     # Reset the springs back to zero for each blow
                elem.u_perm = 0.0
            self.base.u_perm = 0.0

            v0 = hammer.impact_velocity(self.total_mass)
            hammer.reset(v0)

            steps = int(t_per_blow / dt)
            for _ in range(steps):
                self.conservation_of_momentum(dt, hammer=hammer)
                blow_time.append(self.history["time"][-1])
                F_top = self.k_axial * (self.u[1] - self.u[0]) + self.c_axial * (self.v[1] - self.v[0])
                blow_force.append(F_top)

            self.blow_depths.append(self.segment_elevations[-1] + self.u[-1])
            self.blow_histories.append({
                "time": np.array(blow_time) - blow_time[0],
                "force": np.array(blow_force),
            })

            for _ in range(int(rest_time / dt)):
                self.conservation_of_momentum(dt, hammer=None)

    def plot_response(self):
        t = self.history["time"]
        plt.figure(figsize=(12, 6))
        plt.subplot(2, 1, 1)
        plt.plot(t, self.history["top_disp"], label="Top displacement")
        plt.plot(t, self.history["tip_disp"], label="Tip displacement")
        plt.ylabel("Displacement (m)")
        plt.legend()
        plt.grid()

        plt.subplot(2, 1, 2)
        plt.plot(t, self.history["top_vel"], label="Top velocity")
        plt.plot(t, self.history["tip_vel"], label="Tip velocity")
        plt.ylabel("Velocity (m/s)")
        plt.xlabel("Time (s)")
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.show()

    def plot_blow_counts(self):
        if not self.blow_depths:
            print("No blow data to plot.")
            return

        blow_counts = np.arange(1, len(self.blow_depths) + 1)
        depths = np.array(self.blow_depths)

        plt.figure(figsize=(5, 6))
        plt.plot(blow_counts, depths, label="Pile tip depth")
        plt.xlabel("Blow count")
        plt.xlim(0)
        plt.ylabel("Penetration depth (m)")
        # plt.ylim(-30)
        plt.legend()
        plt.show()

    def plot_blow_forces(self, n_to_plot=None):
        plt.figure(figsize=(8, 5))
        if n_to_plot is None:
            n_to_plot = len(self.blow_histories)
        for i in range(n_to_plot):
            data = self.blow_histories[i]
            plt.plot(data["time"] * 1e3, data["force"] * 1e-3, label=f"Blow {i+1}")
        plt.xlabel("Time after blow [ms]")
        plt.ylabel("Force at top [kN]")
        plt.legend()
        plt.title("Force-time at pile top after each blow")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

#$$ === Main Execution ===
if __name__ == "__main__":
    # Define hammer (Junttan HHK3a)
    hammer = HAMMER(mass=3000,   # kg
                    efficiency=0.8,  # energy efficiency
                    E_rated=35000, # J [Nm]
                    cushion_stiffness=5e6, 
                    cushion_damping=0.2, 
                    cushion_quake=0.002)


    # Define soil layers from top down (higher elevation to lower)
    Qs_max1 = 0.005*10*np.pi*(0.5)*7  # alpha_s * qc * area [MN]
    Qs_max2 = 0.005*30*np.pi*(0.5)*7
    layers = [
        SoilLayer(top=0.0, bottom=-7.0, Qs_max=Qs_max1*1e6, damping=0.16, shaft_quake=0.002),     # soft sand
        SoilLayer(top=-7.0, bottom=-15.0, Qs_max=Qs_max2*1e6, damping=0.16, shaft_quake=0.002),    # dense sand
    ]
    profile = SoilProfile(layers)

    # Define pile
    pile = PILE(length=20.0, diameter=0.5, bottom_elevation=-5,
                n_increments=10,
                youngs_modulus=35E3, damping_ratio=0.02, 
                material_density=2500, # kg/m3 (concretw)
                soil_profile=profile)

    # Simulation parameters
    dt = 1e-4
    t_per_blow = 0.05   # simulate 10 ms after impact
    rest_time = 0.5    # 1 second rest between blows
    n_blows = 10

    pile.simulate(hammer=hammer, n_blows=n_blows, dt=dt,
                  t_per_blow=t_per_blow, rest_time=rest_time)

    pile.plot_response()
    pile.plot_blow_counts()
    pile.plot_blow_forces(n_to_plot=None)
