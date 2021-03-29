4 MPCs in this repo : 

Common to all models :  

Optional : 2 boolean (accessible with python binding )     
	-> relative_forces : To add the relatives force in the cost function 
	-> implicit_integration : Implicit integration scheme  
	
The dynamics need to be updated before each control cycle with update_model (with python binding) since      
the dynamics depend on the position of the feet in the predicted time horizon.     

1 - MPC linear :  cost_function.pdf ; cost_shoulder_contact_point.pdf       
-----------------------------------------------------------------------       

The state :   X = [x,y,z,roll,pitch,yaw,vx,vy,vz, v_roll, v_pitch, v_yaw] x12       


--> quadruped Model :    
The command : U = [fx1,fy1,fz1, ...  fx4,fy4,fz4] x12   
(the feet not in contact f = 0)   
 
The cost function consists of : cf cost_function.pdf   

	- Quadratic cost : 0.5*|| X - X_ref  ||^2   
	- Quadratic cost norm of U (relative) : 0.5* || U - 9.81*m/nb_contact ||^2    
	- Shoulder <--> contact point non linear cost : can be set to 0     
	(cf cost_shoulder_contact_point.pdf)    



2 - MPC non linear : non_linear_dynamics.pdf    
------------------------------------------------------------------    

--> quadruped_nl Model :    
Close to MPC linear, the B matrix  and its derivatives are different.    


3 - MPC footstep optimization : footstep_optimization.pdf     
-----------------------------------------------------------    

The state :   X = [x,y,z,roll,pitch,yaw,vx,vy,vz, v_roll, v_pitch, v_yaw    
                   px1,py1, ...   px4,py4] x20  (position of the feet in 2D)    


--> quadruped_augmented Model :    
The command : U = [fx1,fy1,fz1, ...  fx4,fy4,fz4] x12    
Same dynamic, same command, state augmented.     
 
--> quadruped_step :    
The command : U = [Delta px, Delta Py, ...] x4    
For now, a maximum of two feet are modified at the same time.    
Node inserted between the augmented model to modify the position of the feet.    


4 - MPC footstep optimization + time optim : time_optimization.pdf     
----------------------------------------------------------------------   
 
The state :   X = [x,y,z,roll,pitch,yaw,vx,vy,vz, v_roll, v_pitch, v_yaw    
                   px1,py1, ...   px4,py4 , dt] x21     
 (position of the feet in 2D + integration time = period phase / nb of nodes)    

 
--> quadruped_augmented_time Model :    
The command : U = [fx1,fy1,fz1, ...  fx4,fy4,fz4] x12    
Same dynamic, same command, state augmented.     

--> quadruped_step_time :    
The command : U = [Delta px, Delta Py, ...] x4    
Node inserted between the augmented models to modify the position of the feet.    

--> quadruped_time     
The command : U = [dt] x1    
Node inserted between the augmented models to modify the integration time.    







