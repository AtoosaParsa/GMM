# A molecular dynamics simulator for a chain of disks with fixed boundary condition

import numpy as np 
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import math

class switch():
    def evaluate(indices, plot=False):
        ## Experimental Parmeters
        # number of particles
        N=20
        # spring constant for harmonic force law
        K=100
        # area density of particles
        rho=2
        # Mass of particles
        #M=abs(1+randn(1,N)/3)
        N_heavy = 10
        #indices = random.sample(range(0, N), N_heavy)
        massRatio = 10
        #M=np.zeros(N)+1
        M=np.ones(N)
        M[indices] = massRatio
        # Diameter of particles
        Dn=np.ones(N) * math.sqrt(4*massRatio/rho)
        D=np.max(Dn)

        # Drag coefficient
        B=.1  

        # amplitude of shaking
        A=D/200
        # frequency of shaking
        w=2*math.pi*0.95/2

        # size of box
        Lx=N*D
        Ly=D

        # compresssion strain [Lx x]->[Lx x]*(1+ep)
        ep=-2*A

        # Display size
        LyD=Ly

        # total simulation time (short for demo)
        TT=100#256

        ## Physical Parameters
        g=0

        ## Display Parameters
        # plot ?
        plotit=True
        # number of timesteps to skip before plotting
        Nplotskip=50

        ## Simulation Parmeters
        dt=1e-2
        # number of time steps
        Nt=int(TT/dt)

        ## Initial Conditions
        #[x y]=ndgrid(D/2:D:Lx-D/2,D/2:D:Ly-D/2)
        x = np.arange(D/2, Lx-D/2+1, D)
        y = np.arange(D/2, Ly-D/2+1, D)
        #ii=randperm(numel(x),N)
        ii = np.random.permutation(x.size)
        #x=reshape(x(ii),1,N)
        x=x[ii]
        #y=reshape(y(ii),1,N)
        y=np.ones(N)*y

        # Lx compression
        x=x*(1+ep)
        Lx=Lx*(1+ep)

        vx=0*np.random.randn(N)/100
        vy=0*np.random.randn(N)/100
        vx=vx-np.mean(vx)
        vy=vy-np.mean(vy)

        ax_old=0*x
        ay_old=0*y

        x_ini = x
        y_ini = y
        ## Save variables
        # List of quantities to be saved at each time step.
        Ek=np.zeros(Nt)    # Kinetic Energy
        Ep=np.zeros(Nt)    # particle-particle potential
        Ewp=np.zeros((Nt, 4))   # wall-particle potential (1234)=>(LBRT)
        xs=np.zeros((Nt, N))    # x-position
        ys=np.zeros((Nt, N))    # y-position
        vxs=np.zeros((Nt, N))   # x-velocity
        vys=np.zeros((Nt, N))   # y-velocity

        ## Setup Plotting
        # clf;
        # h=zeros(1,N);
        # for np=1:N
        #   %cc=(M(np)-max(M))/(min(M)-max(M))*[1 1 1];
        #   cc=[1 0 0];
        #   h(np)=rectangle('Position',[x(np)-.5*D y(np)-.5*D D D],'Curvature',[1 1],...
        #     'edgecolor','k','facecolor',cc);
        # end
        # hold on;
        # hl=plot([0 0]*Lx,[0 Ly],'k',[1 1]*Lx,[0 Ly],'k',[0 Lx],[0 0]*Ly,'k',[0 Lx],[1 1]*Ly,'k');
        # hold off;
        # axis('equal');
        # axis([-A Lx+A 0 LyD]);
        # %pause;

        ## Main Loop
        Rw=Lx
        Lw=0
        inlw = np.zeros(Nt)
        for nt in range(1, Nt):
  
            # update wall position
            #Rw=Lx+A*sin(w*nt*dt);
            Lw=A*np.sin(w*nt*dt)
            inlw[nt] = Lw
            
            # plot particles
            #   if(plotit && rem(nt-1,Nplotskip)==0)
            #     for np=1:N
            #       set(h(np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)]);
            #     end
            #     set(hl(1),'xdata',Lw*[1 1]);
            #     %set(hl(2),'xdata',Rw*[1 1]);
            #     axis('equal');
            #     axis([-A Lx+A 0 LyD]);
            #     drawnow;
            #   end
            
            # first step in Verlet integration
            x=x+vx*dt+ax_old*dt**2/2
            y=y+vy*dt+ay_old*dt**2/2
            
            # position dependent calculations
            # save positions
            xs[nt,:]=x
            ys[nt,:]=y

            # Interaction detector and Force Law
            Fx=np.zeros(N)
            Fy=np.zeros(N)
            
            for nn in range(1, N):
                for mm in range(nn+1, N):
                    dy=y[mm]-y[nn]
                    Dnm=(Dn[nn]+Dn[mm])/2
                    if(abs(dy) <= Dnm):
                        dx=x[mm]-x[nn]
                        dnm=dx**2+dy**2
                        if(dnm<Dnm**2):
                            dnm=np.sqrt(dnm)
                            F=-K*(Dnm/dnm-1)
                            # particle-particle PE
                            Ep[nt]=Ep[nt]+(Dnm-dnm)**2
                            # particle-particle Force Law
                            Fx[nn]=Fx[nn]+F*dx
                            Fx[mm]=Fx[mm]-F*dx
                            # particle-particle Force Law
                            Fy[nn]=Fy[nn]+F*dy
                            Fy[mm]=Fy[mm]-F*dy
            
            Ep[nt]=K/2*Ep[nt]
            
            ii=x<Lw+Dn/2
            # Left wall
            dw=x[ii]-Lw-Dn[ii]/2
            Fx[ii]=Fx[ii]-K*dw;  
            # PE
            Ewp[nt,0]=K*np.sum(dw**2)/2
            
            ii=y<Dn/2
            # Bottom wall
            dw=y[ii]-Dn[ii]/2
            Fy[ii]=Fy[ii]-K*dw
            # PE
            Ewp[nt,1]=K*np.sum(dw**2)/2
            
            ii=x>Rw-Dn/2
            # Right wall
            dw=x[ii]-(Rw-Dn[ii]/2)
            Fx[ii]=Fx[ii]-K*dw
            # PE
            Ewp[nt,2]=K*np.sum(dw**2)/2
            
            ii=y>Ly-Dn/2
            # Top wall
            dw=y[ii]-(Ly-Dn[ii]/2)
            Fy[ii]=Fy[ii]-K*dw
            # PE
            Ewp[nt,3]=K*np.sum(dw**2)/2
            
            # correction for velocity dependent force
            ax=(Fx/M-B*(vx+ax_old*dt/2))/(1+B*dt/2) 
            ay=(Fy/M-B*(vy+ay_old*dt/2)-g)/(1+B*dt/2)
            
            # second step in Verlet integration
            vx=vx+(ax_old+ax)*dt/2
            vy=vy+(ay_old+ay)*dt/2

            # velocity dependent calculations
            # Kinetic energy
            Ek[nt]=np.sum(M*(vx**2+vy**2))/2
            # save velocities
            vxs[nt,:]=vx
            vys[nt,:]=vy

            ax_old=ax
            ay_old=ay

        fitness=np.sum(np.abs(xs[:, N-1]))
        if (plot):
            plt.figure(figsize=(6.4,4.8))
            plt.plot(list(range(Nplotskip, Nt)), xs[Nplotskip:Nt, N-1], color='blue')
            plt.xlabel("Time Steps")
            plt.ylabel("Position")
            plt.title("Position of the Last Particle", fontsize='small')
            plt.grid(color='skyblue', linestyle=':', linewidth=0.5)
            plt.tight_layout()
            plt.show()

            plt.figure(figsize=(6.4,4.8))
            plt.plot(list(range(Nplotskip, Nt)), inlw[Nplotskip:Nt], color='blue')
            plt.xlabel("Time Steps")
            plt.ylabel("Position")
            plt.title("Position of the First Particle", fontsize='small')
            plt.grid(color='skyblue', linestyle=':', linewidth=0.5)
            plt.tight_layout()
            plt.show()

        return fitness

    def showPacking(indices, save=0):
        ## Experimental Parmeters
        # number of particles
        N=20
        # spring constant for harmonic force law
        K=100
        # area density of particles
        rho=2
        # Mass of particles
        #M=abs(1+randn(1,N)/3)
        N_heavy = 10
        #indices = random.sample(range(0, N), N_heavy)
        massRatio = 10
        #M=np.zeros(N)+1
        M=np.ones(N)
        M[indices] = massRatio
        # Diameter of particles
        Dn=np.sqrt(4*M/rho)
        D=np.max(Dn)

        # size of box
        Lx=N*D
        Ly=D

        # Drag coefficient
        B=.1  

        # amplitude of shaking
        A=D/200
        # frequency of shaking
        w=2*math.pi*0.95/2

        # size of box
        Lx=N*D
        Ly=D

        # compresssion strain [Lx x]->[Lx x]*(1+ep)
        ep=-2*A

        # Display size
        LyD=Ly

        # total simulation time (short for demo)
        TT=256

        ## Physical Parameters
        g=0

        ## Display Parameters
        # plot ?
        plotit=True
        # number of timesteps to skip before plotting
        Nplotskip=50

        ## Simulation Parmeters
        dt=1e-2
        # number of time steps
        Nt=int(TT/dt)

        ## Initial Conditions
        #[x y]=ndgrid(D/2:D:Lx-D/2,D/2:D:Ly-D/2)
        x = np.arange(D/2, Lx-D/2+1, D)
        y = np.arange(D/2, Ly-D/2+1, D)
        #ii=randperm(numel(x),N)
        ii = np.random.permutation(x.size)
        #x=reshape(x(ii),1,N)
        x=x[ii]
        #y=reshape(y(ii),1,N)
        y=np.ones(N)*y

        # Lx compression
        x=x*(1+ep)
        Lx=Lx*(1+ep)

        x_ini = x
        y_ini = y

        m_min = min(M)
        m_max = max(M)
        fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
        ells = []
        m_all = []
        for i in range(N):
            x_now = x_ini[i]%Lx
            y_now = y_ini[i]%Ly
            for k in range(-1, 2):
                for l in range(-1, 2):                        
                    e = Ellipse((x_now+k*Lx, y_now+l*Ly), D,D,0)
                    ells.append(e)
                    m_all.append(M[i])

        i = 0
        for e in ells:
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_facecolor('C2')
            
            e.set_alpha(0.2+(m_all[i]-m_min)/(m_max-m_min+1)*0.8)
            e.set_alpha(0.2+(m_all[i]-m_min)/(m_max-m_min+1)*0.3)

            i += 1
                    
        ax.set_xlim(0, Lx)
        ax.set_ylim(0, Ly)

        plt.show() 
        if save == 1:
            fig.savefig(hn, dpi = 300)