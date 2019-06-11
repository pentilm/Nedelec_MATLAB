function Solve( obj )
% Organize the whole algorithm
Initialize( obj );
Partition2D( obj );
PreCal( obj );
Form_Mat( obj );
Init_Val( obj );
Marching( obj );
Post_Process( obj );
% Plot_Sol( obj );
end

