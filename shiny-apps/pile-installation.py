# Full Python script with two piles, animation, and updated resistance tracking

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.gridspec import GridSpec
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch



# Parameters
num_blows = 10
initial_velocity = 5.0       # Initial velocity imparted by hammer [m/s]
blow_duration = 0.3
frames_per_blow = 30
frame_dt = blow_duration / frames_per_blow

# Calibration chamber parameters
soil_depth = 20.0
layer_depths = [14.0, 20.0]
starting_penetration = 4.0

pile_length = 15.0
num_nodes = 10
monitor_node_index = num_nodes - 2  # near the tip

# Pile A parameters
params_A = {
    'shaft_stiffness': [80, 120],
    'base_resistance': 400,
    'damping': 15.0,
    'x_offset': -0.5,
    'color': 'blue',
    'label': 'Pile A'
}

# Pile B parameters
params_B = {
    'shaft_stiffness': [80, 120],
    'base_resistance': 400,
    'damping': 30.0,
    'x_offset': 0.5,
    'color': 'green',
    'label': 'Pile B'
}

# Pile class
class PILE:
    def __init__(self, length, tip_depth, num_nodes, shaft_stiffness, base_resistance, damping, x_offset):
        self.length = length
        self.tip_depth = tip_depth
        self.num_nodes = num_nodes
        self.shaft_stiffness = shaft_stiffness
        self.base_resistance = base_resistance
        self.damping = damping
        self.x_offset = x_offset
        self.update_node_positions()
        
        # Set up soil
        self.cpt = pd.DataFrame()    # z=0m = seabed. Negative downwards
        self.layers = {}

    def update_node_positions(self):
        self.nodes_y = np.linspace(self.tip_depth - self.length, self.tip_depth, self.num_nodes)
        self.nodes_x = np.full(self.num_nodes, self.x_offset)

    def get_nodes(self):
        return list(zip(self.nodes_x, self.nodes_y))

    def get_spring_constants(self):
        segment_length = self.length / (self.num_nodes - 1)
        spring_constants = np.zeros(self.num_nodes - 1)
        for i in range(self.num_nodes - 1):
            node_mid_depth = (self.nodes_y[i] + self.nodes_y[i + 1]) / 2
            if node_mid_depth <= layer_depths[0]:
                k = self.shaft_stiffness[0] * segment_length
            else:
                k = self.shaft_stiffness[1] * segment_length
            spring_constants[i] = k
        return spring_constants

    def _static_resistance(self, current_depth):
        # Shaft resistance
        Fs_all_layers = []
        for layer in self.layers:       # Cycle through each layer along the pile
            z_top_inc = self.layers[layer][0]
            z_base_inc = self.layers[layer][-1]
            if z_top_inc < current_depth:           # If layer is below pile, ignore.
                continue
            
            else:
                current_layer_id = layer 
                
                if z_base_inc < current_depth:
                    z_base_inc = current_depth
                z_mid_inc = (z_base_inc + z_top_inc)/2
                    
                # Assuming just sand only
                qc_mean = self.cpt.loc[(self.cpt.z > z_base_inc) & (self.cpt.z > z_top_inc)].qc.mean()
                Fs = 0.005 * qc_mean * np.pi * self.D * (z_top_inc - z_base_inc)    #[MN]
                Fs_all_layers.append(Fs)
                
        Fs_total = sum(Fs_all_layers)
            
        # Base resistnace
        near_z_ix = self.cpt.z.sub(current_depth).abs().idxmin()    # Index of row with depth closest to pile depth
        qt_tip = self.cpt.qt.iloc[near_z_ix]                    # qt at pile tip [MPa]
        alpha_p = 0.5
        Fbear = alpha_p * qt_tip * ((np.pi*self.D**2)/4)     # MN

        return Fs_total + Fbear    # [MN]

    def _dynamic_resistance(self, current_depth, velocity):

    def _system_efficiency(self, HAMMER):
        e = 0.5         # Coefficient of resistution (for steel)
        m_ram = HAMMER.ram_mass
        mp = self.mass
        if m_ram >= e*mp:
            nu = (m_ram + (e**2)*mp)/(m_ram + mp)
        elif m_ram < e*mp:
            nu = (m_ram + (e**2)*mp)/(m_ram + mp) - ((m_ram + e*mp)/(m_ram + mp))**2
        
        return nu

    def initial_velocity(self,HAMMER):
        """Velocity created by kinetic energy of hammer blow"""
        nu = self._system_efficiency(HAMMER)
        nu = 1      # Says system efficiency is equal to 0.02??? Does this make sense?
        
        self.received_blows+=1
        E_rated = HAMMER.E_rated*1E3  # KJ to J
        
        v0 = np.sqrt((2*nu*HAMMER.efficiency*E_rated)/(self.mass + HAMMER.mass))
        
        return v0*-1       # [m/s]
    
    
    
    def simulate_pile_movement(self, HAMMER, initial_velocity, max_depth, time_step):
        tolerance = 0.001  # If velocity in this range: then pile has stopped [m/s]
        
        current_depth = self.ztip
        velocity = initial_velocity
    
        Wh = (HAMMER.mass * 9.81)/1E6  # Weight of hammer [MN]
        Wh = self.hammer_on_pile * Wh   # Is hammer (partially) on pile? Then account for this 
        Wf = (HAMMER.follower_mass * 9.81)/1E6 # Weight of follower [MN]
        Wf = self.follower_on_pile * Wf
        mass_on_pile = self.hammer_on_pile * HAMMER.mass + self.follower_on_pile * HAMMER.follower_mass     # [kg]
        Wp = (self.mass * 9.81)/1E6   # Weight of pile [MN]
    
        iteration = 0
        while current_depth > max_depth:
            # print(iteration)
            SRD = self._static_resistance(current_depth,velocity)
            
            Fbuoy = self._soil_buoyancy(current_depth) + self._buoyancy(current_depth)
            Fdrag = self._inertial_drag(HAMMER,velocity)    

            net_force = (Fbuoy + Fshaft + Fbear + Fdrag) - (Wp + Wh + Wf)       # [MN]. Negative equilibrium sends pile downwards
            srd = Fbear + Fshaft    # To be saved for the results
            acceleration = (net_force*1E6) / (self.mass + mass_on_pile)   # [m/s2]
            
            updated_velocity = velocity + (acceleration * time_step)
            
            if updated_velocity > tolerance:  # Pile can't move upwards. Reduce time_step & try again
                time_step*=0.1   
                continue
                
            else:       # Update everything and check stopping condition
                velocity += acceleration * time_step
                current_depth += velocity * time_step
                self.ztip = current_depth
                iteration += 1
                
                rate_factor = np.nan # Cant plot this now that it's variable by soil layer. DUD for now
                self.results.loc[len(self.results)] = [current_depth,Fbear,Fshaft,Fbear_unfactored,Fshaft_unfactored,
                                                       Fbuoy,Fdrag,srd,velocity*-1,
                                                       self.received_blows,rate_factor,norm_velocity_ann,norm_velocity_Douter]    # Append results. Blows = 0 
                
                # Check for stopping condition
                if self.stopping_condition_met(current_depth, velocity, max_depth,tolerance):
                    break
    
        # return
    

# Physics functions
def velocity_t(t, v0, damping):
    return v0 * np.exp(-damping * t)

def penetration_increment(v0, t, damping):
    return (v0 / damping) * (1 - np.exp(-damping * t))

def resistance_at_time(pile, penetration_start, t, v0, node_index):
    penetration_now = penetration_start + penetration_increment(v0, t, pile.damping)
    spring_constants = pile.get_spring_constants()
    node_top = pile.nodes_y[node_index]
    node_bottom = pile.nodes_y[node_index + 1]
    node_mid = (node_top + node_bottom) / 2
    k = spring_constants[node_index]
    embed_length = max(0, penetration_now - node_mid)
    velocity = velocity_t(t, v0, pile.damping)
    resistance = k * embed_length + pile.damping * velocity
    return resistance, penetration_now, velocity

# Filter only parameters expected by Pile
physical_keys = ['shaft_stiffness', 'base_resistance', 'damping', 'x_offset']
pileA_args = {k: params_A[k] for k in physical_keys}
pileB_args = {k: params_B[k] for k in physical_keys}

pileA = PILE(pile_length, starting_penetration, num_nodes, **pileA_args)
pileB = PILE(pile_length, starting_penetration, num_nodes, **pileB_args)

# Plot setup
fig = plt.figure(figsize=(14, 8))
gs = GridSpec(3, 2, width_ratios=[1, 1], height_ratios=[1, 1, 1], figure=fig, wspace=0.3, hspace=0.4)

ax0 = fig.add_subplot(gs[:, 0])
ax0.set_xlim(-1.5, 1.5)
ax0.set_ylim(soil_depth + 1, -1)
ax0.set_title('Pile Penetration Animation')
ax0.fill_between([-2, 2], 0, layer_depths[0], color='yellow', alpha=0.3)
ax0.fill_between([-2, 2], layer_depths[0], layer_depths[1], color='orange', alpha=0.3)

def init_pile_plot(pile, color):
    nodes = pile.get_nodes()
    nodes_x, nodes_y = zip(*nodes)
    pile_line, = ax0.plot(nodes_x, nodes_y, color=color, linewidth=6)
    nodes_plot, = ax0.plot(nodes_x, nodes_y, 'ko', markersize=6)
    springs_lines = []
    for i in range(len(nodes_x) - 1):
        spring_line, = ax0.plot([nodes_x[i], nodes_x[i+1]], [nodes_y[i], nodes_y[i+1]], 'gray', linestyle='--')
        springs_lines.append(spring_line)
    arrow = FancyArrowPatch(
        (pile.x_offset, pile.tip_depth - pile_length / 2),
        (pile.x_offset, pile.tip_depth - pile_length / 2 - 0.5),
        arrowstyle='-|>', mutation_scale=20, color='red', linewidth=2
    )
    ax0.add_patch(arrow)
    return pile_line, nodes_plot, springs_lines, arrow

pileA_artists = init_pile_plot(pileA, params_A['color'])
pileB_artists = init_pile_plot(pileB, params_B['color'])

# Resistance vs time
ax1 = fig.add_subplot(gs[0:2, 1])
ax1.set_xlim(0, blow_duration)
ax1.set_ylim(0, 1000)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Resistance (kN)')
ax1.set_title('Resistance vs Time at Node')
resA_line, = ax1.plot([], [], '-', color=params_A['color'], label='Pile A')
resB_line, = ax1.plot([], [], '-', color=params_B['color'], label='Pile B')
resA_point, = ax1.plot([], [], 'o', color=params_A['color'])
resB_point, = ax1.plot([], [], 'o', color=params_B['color'])
node_label = ax1.text(0.01, 0.95, '', transform=ax1.transAxes, verticalalignment='top')
ax1.legend()

# Resistance vs depth
ax2 = fig.add_subplot(gs[2, 1])
depth_samples = np.linspace(0, soil_depth, 200)
shaftA = (np.minimum(depth_samples, layer_depths[0]) * params_A['shaft_stiffness'][0] +
          np.maximum(depth_samples - layer_depths[0], 0) * params_A['shaft_stiffness'][1])
shaftB = (np.minimum(depth_samples, layer_depths[0]) * params_B['shaft_stiffness'][0] +
          np.maximum(depth_samples - layer_depths[0], 0) * params_B['shaft_stiffness'][1])
baseA = np.where(depth_samples >= pile_length, params_A['base_resistance'], 0)
baseB = np.where(depth_samples >= pile_length, params_B['base_resistance'], 0)
ax2.plot(shaftA + baseA, depth_samples, label='Pile A', color=params_A['color'])
ax2.plot(shaftB + baseB, depth_samples, label='Pile B', color=params_B['color'])
tipA_line = ax2.axhline(pileA.tip_depth, color=params_A['color'], linestyle='--')
tipB_line = ax2.axhline(pileB.tip_depth, color=params_B['color'], linestyle='--')
ax2.invert_yaxis()
ax2.set_xlabel('Resistance (kN)')
ax2.set_ylabel('Depth (m)')
ax2.set_title('Static Resistance vs Depth')
ax2.legend()

# Animation update
resA_hist = []
resB_hist = []
penetration_A = [starting_penetration]
penetration_B = [starting_penetration]
time_vec = np.linspace(0, blow_duration, frames_per_blow)

def update(frame):
    blow_idx = frame // frames_per_blow
    frame_in_blow = frame % frames_per_blow
    t_now = time_vec[frame_in_blow]
    if blow_idx >= num_blows:
        return []

    pen_start_A = penetration_A[blow_idx]
    pen_start_B = penetration_B[blow_idx]
    pen_now_A = pen_start_A + penetration_increment(initial_velocity, t_now, pileA.damping)
    pen_now_B = pen_start_B + penetration_increment(initial_velocity, t_now, pileB.damping)
    pileA.tip_depth = pen_now_A
    pileB.tip_depth = pen_now_B
    pileA.update_node_positions()
    pileB.update_node_positions()
    if frame_in_blow == frames_per_blow - 1:
        penetration_A.append(pen_now_A)
        penetration_B.append(pen_now_B)

    for pile, artists in zip([pileA, pileB], [pileA_artists, pileB_artists]):
        pile_line, nodes_plot, springs_lines, arrow = artists
        nodes = pile.get_nodes()
        nodes_x, nodes_y = zip(*nodes)
        pile_line.set_data(nodes_x, nodes_y)
        nodes_plot.set_data(nodes_x, nodes_y)
        for i, spring_line in enumerate(springs_lines):
            spring_line.set_data([nodes_x[i], nodes_x[i+1]], [nodes_y[i], nodes_y[i+1]])
        arrow.set_alpha(1 if frame_in_blow < 5 else 0)  # only show arrow briefly
        arrow.set_positions(
            (pile.x_offset, pile.tip_depth - pile_length / 2),
            (pile.x_offset, pile.tip_depth - pile_length / 2 - 0.5)
        )
        arrow.set_visible(frame_in_blow < 5)

    rA, _, _ = resistance_at_time(pileA, pen_start_A, t_now, initial_velocity, monitor_node_index)
    rB, _, _ = resistance_at_time(pileB, pen_start_B, t_now, initial_velocity, monitor_node_index)
    if frame_in_blow == 0:
        resA_hist.clear()
        resB_hist.clear()
    resA_hist.append((t_now, rA))
    resB_hist.append((t_now, rB))
    timesA, valsA = zip(*resA_hist)
    timesB, valsB = zip(*resB_hist)
    resA_line.set_data(timesA, valsA)
    resB_line.set_data(timesB, valsB)
    resA_point.set_data(t_now, rA)
    resB_point.set_data(t_now, rB)
    tipA_line.set_ydata(pileA.tip_depth)
    tipB_line.set_ydata(pileB.tip_depth)
    node_depth = (pileA.nodes_y[monitor_node_index] + pileA.nodes_y[monitor_node_index + 1]) / 2
    node_label.set_text(f'Monitor Node Depth: {node_depth:.2f} m')
    return []

# Run animation
anim = FuncAnimation(fig, update, frames=num_blows * frames_per_blow, interval=50, blit=False)
plt.show()
