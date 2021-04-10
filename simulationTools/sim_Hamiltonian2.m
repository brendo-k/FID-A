%sim_Hamiltonian.m
%Robin Simpson and Jamie Near, 2014.
%Kimberly Chan added separate Hamiltonian for J-coupling only, for use
%during shaped rf pulses (sim_shapedRF.m).
%
% USAGE:
% [H,d] = sim_Hamiltonian(sys,Bfield);
%
% DESCRIPTION:
% Creates the nxn Hamiltonian matrix for a spin system which can then be
% used in other functions to simulate NMR experiments.
%
% INPUTS:
% sys     = spin system definition structure.
% Bfield  = magnetic field strength (Tesla).
%
% OUTPUTS:
% H       = n x n Hamiltonian matrix for spin system.
% d       = Equilibrium density matrix.

function [H,d] = sim_Hamiltonian2(sys,Bfield)

    omega0 = -2*pi*Bfield*42576000;

    for n=1:length(sys) %JN - Looping through the parts of the spin system:
        %initialize parameters:
        nspins = length(sys(n).shifts);

        %The J matrix and the shifts vector have to be
        %converted into radians and rad/s, respectively:
        H(n).J=sys(n).J*2*pi;
        H(n).shifts=sys(n).shifts;
        H(n).shifts_rads = sys(n).shifts*(omega0/10^6);

        %add some other "header" information to the Hamiltonian structure;
        H(n).nspins=nspins;
        H(n).Bfield=Bfield;

        H(n).Fz = zeros(2^nspins);
        H(n).Fy = zeros(2^nspins);
        H(n).Fx = zeros(2^nspins);
        H(n).HAB = zeros(2^nspins,2^nspins);
        H(n).Ix = zeros(2^nspins,2^nspins,nspins);
        H(n).Iy = zeros(2^nspins,2^nspins,nspins);
        H(n).Iz = zeros(2^nspins,2^nspins,nspins);

        d{n}=H(n).Fz;

        %Store the spin-system scaling factor in the H structure, and scale the
        %density matrix for each part of the spin system:
        H(n).scaleFactor=sys(n).scaleFactor;
        d{n}=d{n}*sys(n).scaleFactor;

        %Basis states
        I0=[1 0;0 1];
        Ix=0.5*[0 1;1 0];
        Iy=(1i/2)*[0 1;-1 0];
        Iz=(1/2)*[-1 0;0 1];

        for spin = 1:nspins
            if spin == 1
                temp_Ix = Ix;
                temp_Iy = Iy;
                temp_Iz = Iz;
            else
                temp_Ix = I0;
                temp_Iy = I0;
                temp_Iz = I0;
            end
            for pos = 2:nspins
                if spin == pos
                    temp_Ix = kron(temp_Ix, Ix);
                    temp_Iy = kron(temp_Iy, Iy);
                    temp_Iz = kron(temp_Iz, Iz);
                else
                    temp_Ix = kron(temp_Ix, I0);
                    temp_Iy = kron(temp_Iy, I0);
                    temp_Iz = kron(temp_Iz, I0);
                end
            end
            H(n).Ix(:,:,spin) = temp_Ix;
            H(n).Iy(:,:,spin) = temp_Iy;
            H(n).Iz(:,:,spin) = temp_Iz;
        end
        H(n).Fx = sum(H(n).Ix, 3);
        H(n).Fy = sum(H(n).Iy, 3);
        H(n).Fz = sum(H(n).Iz, 3);

            %Now the resonance component - each spin has a resonance frequency
            %omega0-chemshift. This component is resfreq*H.Iz for each spin. In fact
            %we assume we are in a frame rotating at omega0 so we only have to
            %include H.shifts*H.Iz for each spin. This means we are not drastically
            %undersampled as we would otherwise be.
            %     for e=1:nspins
            %         rescomp = (H(n).shifts_rads(e))*H(n).basis(:,:,e);
            %         H(n).HAB(:,:) = H(n).HAB(:,:) + diagonal(:,:).*rescomp;
            %     end
    end
end




