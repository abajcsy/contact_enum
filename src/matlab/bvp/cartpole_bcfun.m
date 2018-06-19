function bc = cartpole_bcfun(ddsdswitch_a,ddthetadswitch_a,dp1dswitch_b,dp2dswitch_b,dp3dswitch_b,dp4dswitch_b,ds0,ds_a,ds_b,dsdswitch_a,dtheta0,dtheta_a,dtheta_b,dthetadswitch_a,l,lmda1_a,lmda2_a,lmda3_a,p1_b,p2_b,p3_b,p4_b,s0,s_a,s_b,sgoal,swall,theta0,theta_a,theta_b,wg)
%CARTPOLE_BCFUN
%    BC = CARTPOLE_BCFUN(DDSDSWITCH_A,DDTHETADSWITCH_A,DP1DSWITCH_B,DP2DSWITCH_B,DP3DSWITCH_B,DP4DSWITCH_B,DS0,DS_A,DS_B,DSDSWITCH_A,DTHETA0,DTHETA_A,DTHETA_B,DTHETADSWITCH_A,L,LMDA1_A,LMDA2_A,LMDA3_A,P1_B,P2_B,P3_B,P4_B,S0,S_A,S_B,SGOAL,SWALL,THETA0,THETA_A,THETA_B,WG)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    18-Jun-2018 16:25:26

t2 = s_b.*2.0;
t3 = sgoal.*2.0;
t4 = t2-t3;
t5 = cos(theta_b);
t6 = l.*lmda1_a.*t5;
bc = [-s0+s_a;-theta0+theta_a;-ds0+ds_a;-dtheta0+dtheta_a;lmda1_a+p1_b-t4.*wg;p2_b+t6;-lmda2_a+p3_b;-lmda3_a+p4_b;dsdswitch_a;dthetadswitch_a;ddsdswitch_a;ddthetadswitch_a;dp1dswitch_b+lmda1_a-t4.*wg;dp2dswitch_b+t6;dp3dswitch_b-lmda2_a;dp4dswitch_b-lmda3_a;-s_b+swall-l.*sin(theta_b);ds_b;dtheta_b];
