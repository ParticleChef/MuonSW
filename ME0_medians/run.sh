executable = run.sh 
universe   = vanilla 
log        = log.log 
output     = out.out 
error      = err.err 
transfer_input_files = x_file.C 
should_transfer_files = YES 
when_to_transfer_output = ON_EXIT
requirements = ( HasSingularity == true)
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el6:latest"
+SingularityBind = "/cvmfs, /cms, /share"
queue


