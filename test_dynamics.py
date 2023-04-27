from math import *
import numpy as np
from icecream import ic
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.pyplot as plt
import matplotlib
import time

def tt_dynamics(params, state, u):
    
    tractor_wheelbase, tractor_hitch_offset, trailer_wheelbase = params[0], params[1], params[2]
    
    # Variable to store the rate of change of state
    alpha = state[4]
    gear = 1
    velocity = state[5]

    if state[3]>pi:
        beta = -1*(2*pi - state[3])
    else:
        beta = state[3]
    
    if(alpha!=0):
        r1 = tractor_wheelbase/np.tan(alpha)
        psi = -1*(alpha/np.abs(alpha))*np.arctan(tractor_hitch_offset/np.abs(r1))
        # print(psi)
        # print(velocity)

        if velocity == 0:
            va = 0
        else:
            va = (velocity/np.abs(velocity))*np.abs(((gear*velocity*tractor_hitch_offset*np.tan(alpha))/(tractor_wheelbase*np.sin(psi))))
        vb = va*np.abs(np.cos(psi-beta))
        r2 = tractor_hitch_offset/np.sin(psi)
        r3 = trailer_wheelbase/np.sin(psi-beta)
        omega_1 = gear*velocity*(1/r1)
        omega_2 = gear*velocity*((r2)/(r1*r3))
        trailer_theta_dot = omega_2
        beta_dot = gear*velocity*(((r2)/(r1*r3)) - (1/r1))

    # Special case when alpha equals zero
    else:
        va = gear*velocity
        psi = 0
        phi = (pi/2) - (np.abs(beta) - (pi/2))
        if(np.sin(phi)==0):
            beta_dot = 0
            trailer_theta_dot = 0
            vb = va*np.abs(np.cos(psi-beta))
        else:
            hitch_turn_radius_zero_alpha = trailer_wheelbase/np.sin(phi)
            trailer_theta_dot = va/hitch_turn_radius_zero_alpha
            beta_dot = trailer_theta_dot
            vb = va*np.abs(np.cos(psi-beta))
    
    q_dot = np.array([vb*np.cos(state[2]), vb*np.sin(state[2]), trailer_theta_dot, beta_dot, u[0], u[1]])
    
    return q_dot


def rk4_dyn(params, x, u):
    
    # RK4 integration with zero-order hold on u
    h = params[3]
    f1 = tt_dynamics(params, x, u)
    f2 = tt_dynamics(params, x + 0.5*h*f1, u)
    f3 = tt_dynamics(params, x + 0.5*h*f2, u)
    f4 = tt_dynamics(params, x + h*f3, u)
    
    next_state = x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
    return next_state

def make_frame(i, params, tractor_coords, trailer_coords, ax):

    # if not includePrevious:
    plt.clf()

    tractor_line_x = np.array([tractor_coords[i,0], tractor_coords[i,0] + (params[0] + params[1])*np.cos(tractor_coords[i,2])])
    tractor_line_y = np.array([tractor_coords[i,1], tractor_coords[i,1] + (params[0] + params[1])*np.sin(tractor_coords[i,2])])


    trailer_line_x = np.array([trailer_coords[i,0], trailer_coords[i,0] + params[2]*np.cos(trailer_coords[i,2])])
    trailer_line_y = np.array([trailer_coords[i,1], trailer_coords[i,1] + params[2]*np.sin(trailer_coords[i,2])])

    # wheel_x = 
    # wheel_y = 
    # wheel_theta = 

    plt.xlim(-20,50)
    plt.ylim(-20,50)
    plt.gca().set_aspect('equal')
    # plt.axis('off')
    artists = []
    artists.append(plt.plot(trailer_line_x, trailer_line_y, color = 'blue', linewidth=5))
    artists.append(plt.plot(tractor_line_x, tractor_line_y, color = 'red', linewidth=5))
    artists.append
    # artists.append(plt.legend(loc='upper right'))   
    print("i value: ", i)

    # if i%10==0:
    #     print("Tractor Coords: ", tractor_coords[i,:])
    #     print("Tractor X: ", tractor_line_x)
    #     print("Tractor Y: ", tractor_line_y)
    #     print("Trailer Coords: ", trailer_coords[i,:])
        # time.sleep(2)


    return artists

def make_gif(params, trajectory, fps, output_path):

    # Extract tractor coords
    tractor_coords = trajectory[:,6:9]
    # Extract trailer coords
    trailer_coords = trajectory[:,0:3]


    fig,ax = plt.subplots()
    numFrames = trajectory.shape[0]
    ani = FuncAnimation(fig, make_frame, repeat=False, frames=numFrames, fargs=[params, tractor_coords, trailer_coords, ax])    
    ani.save(output_path, dpi=50, writer=PillowWriter(fps=fps))
    print("Saved gif to: ", output_path)

def expand_state(params, state):

    tractor_theta = 0
    tractor_wheelbase, tractor_hitch_offset, trailer_wheelbase = params[0], params[1], params[2]
    
    # if state[3]<=pi:
    #     tractor_theta = pi - state[3] + (pi/2) #state[3] + (pi - state[4]);
    # else:
    #     tractor_theta = (pi/2) - (pi - state[3]) #+ (pi/2) #state[3] + (pi - state[4]);

    tractor_theta = pi - state[3] + state[2]
    
    # if state[3]>2*pi:
    #     print("NEGATIVE")

#     tractor_theta = (state[3] - (pi - abs(state[4]))); 
        
    tractor_axle_x = state[0] + trailer_wheelbase*np.cos(state[2]) + tractor_hitch_offset*np.cos(tractor_theta)
    tractor_axle_y = state[1] + trailer_wheelbase*np.sin(state[2]) + tractor_hitch_offset*np.sin(tractor_theta)
    tractor_hitch_x = state[0] + trailer_wheelbase*np.cos(state[2])
    tractor_hitch_y = state[1] + trailer_wheelbase*np.sin(state[2])
        
    return tractor_hitch_x, tractor_hitch_y, tractor_theta




if __name__ == "__main__":

    #tractor wheelbase, trailer wheelbase, tractor hitch offset, dt
    tractor_wheelbase = 3
    tractor_hitch_offset = 0.5
    trailer_wheelbase = 6
    dt = 0.05
    params = (tractor_wheelbase, tractor_hitch_offset, trailer_wheelbase, dt)

    q_init = np.array([0,0,pi/2, pi, 0, 0])
    u = np.array([0.1, 0.1])

    q_dot = tt_dynamics(params, q_init, u)



    # Read trajectory from a text file



    # Get trajectory from arbitrary controls
    N = 500
    init_state = np.array([1,2,pi/2,pi,0,0])
    controls = np.zeros((N-1, 2))
    # controls[50:200,0] = 0.05
    # controls[200:400, 0] = -0.05
    controls[:100,1] = 0.3
    controls[350:,1] = -0.1

    # Make the controls
    max_steer = 0.45
    min_steer = -0.45
    ctrl = -0.05
    steer_sum = 0
    for k in range(N-1):
        controls[k,0] = ctrl

        steer_sum = steer_sum + controls[k,0]*dt
        # print(steer_sum)
        if steer_sum<min_steer:
            ctrl = 0

    # print(controls[:,0])
    # print(np.count_nonzero(controls[:,0]))    
    
    test_trajectory = np.zeros((N,9))
    test_trajectory[0,:6] = init_state
    for j in range(1,N):
        test_trajectory[j,:6] = rk4_dyn(params, test_trajectory[j-1,:6], controls[j-1,:])
        # print(isnan(test_trajectory[j,0]))
        # if isnan(test_trajectory[j,0]):
        #     break
    np.savetxt("./output_traj_test1", test_trajectory, fmt='%.18e', delimiter=' ', newline='\n')
    
    # exit()

    # Expand all the states in the trajectory
    for k in range(N):
        x, y, theta = expand_state(params, test_trajectory[k,:])
        test_trajectory[k,6:9] = [x, y, theta]

    # Random trajectory
    # trajectory = np.array([[5,3,pi/2, pi,0,0,5,9,pi/2],[5,4,pi/2, pi,0,0,5,10,pi/2],[5,5,pi/2, pi,0,0,5,11,pi/2]])
    
    # print(test_trajectory[0,:])
    
    fps = 1/dt
    output_path = "./test_gif.gif"
    make_gif(params, test_trajectory, fps, output_path)

    print(q_dot)
