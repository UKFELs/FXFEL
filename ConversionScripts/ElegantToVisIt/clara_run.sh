dirname=$1_$2
run=1
png=0

if [ ! -d $dirname ];then
    mkdir $dirname
fi

cp clara_run.sh $dirname/clara_run.sh

if [ $run -eq 1 ];then
    cp $2 $dirname/$2
    cd $dirname
    astra2elegant $2 $2.SDDS
    sddsanalyzebeam $2.SDDS $2.analyse.SDDS
    for INP in ../Sz*
      do
      newname=`basename $INP`
      cp ../$newname $newname >& /dev/null
    done
    for INP in ../Sx*
      do
      newname=`basename $INP`
      cp ../$newname $newname >& /dev/null
    done
    ln -s ../count_xyz.tcl count_xyz.tcl >& /dev/null
    cp ../clara_V8_rfcw.lte clara_V8_rfcw.lte
    cp ../clara.track.ele clara.track.ele
    time elegant clara.track.ele > clara.track_out
    tclsh count_xyz.tcl
    cd ../
fi

if [ -e $dirname/clara.ran ];then 
    echo "Performing analysis in directory $dirname - normalising to 10k particles" 
    
    sddsplot -device=cpostscript -layout=1,2 -output=$dirname/clara_cx.ps -lspace=0.45,0.52,-0.51,0.49 -linetypedefault=0,thickness=2 -title= -pSpace=0.15,0.93,0.15,0.93 -graph=line,vary -split=pages $dirname/clara.cen -column=s,\(Cx,Cy\) -column=s,Profile -endpanel -overlay=xmode=norm,yfact=0.04,yoffset=0 -graph=line,type=0 $dirname/clara.mag >& /dev/null
    sddsplot -device=cpostscript -output=$dirname/clara_lattice.ps -lspace=0.45,0.52,-0.51,0.79 -layout=1,2 -linetypedefault=0,thickness=2 -title= -pSpace=0.15,0.93,0.15,0.93 -graph=line,vary -split=pages -toPage=1 $dirname/clara.twi -column=s,beta? -scales=0,0,-10,60 -legend -column=s,etax -yscalesGroup=id=eta -legend -column=s,Profile -endpanel -overlay=xmode=norm,yfact=0.04,yoffset=0 -graph=line,type=0 $dirname/clara.mag -scales=0,0,-10,60 -column=s,\(Sx,Sy\) -graph=line,vary -factor=yMultiplier=1e3 -ylabel="\$gs\$r (mm)" $dirname/clara.sig -split=pages -toPage=3 -column=s,Profile -overlay=xmode=norm,yfact=0.04,yoffset=0 -graph=line,type=0 $dirname/clara.mag -scales=0,0,-0.1,1.5 >& /dev/null
	
    sddsplot -device=cpostscript -output=$dirname/clara_floor.ps -lspace=0.45,0.52,-0.51,0.49 -layout=1,2 -linetypedefault=0,thickness=2 -title= -pSpace=0.15,0.93,0.15,0.93 -graph=line,vary $dirname/clara.floor -column=Z,X -endpanel -scales=0,0,-1,1 -column=s,pCentral $dirname/clara.cen -split=pages -factor=yMultiplier=0.511 -ylabel="p (MeV/c)" -column=s,Profile -overlay=xmode=norm,yfact=0.04,yoffset=0 -mode=y=linear -graph=line,type=0 $dirname/clara.mag -scales=0,0,-20,300 >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/clara_r56.ps -lspace=0.15,0.22,-1.0,0.5 -layout=1,2 -linetypedefault=0,thickness=2 -title= -pSpace=0.15,0.93,0.15,0.93 -graph=line,type=3 -split=pages $dirname/clara.matr -column=s,R56 -legend -column=s,T566 -graph=line,vary=type=2 -legend -column=s,Profile -endpanel -overlay=xmode=norm,yfact=0.04,yoffset=0 -graph=line,type=0 $dirname/clara.mag -scales=0,0,-0.1,0.03 -column=s,Ss -factor=yMultiplier=1e3 $dirname/clara.sig -ylabel="\$gs\$r\$bz\$n (mm)" -column=s,Profile -overlay=xmode=norm,yfact=0.04,yoffset=1e-4 -graph=line,type=0 $dirname/clara.mag -scales=0,0,-0.05,0.4 >& /dev/null

    sddsplot -device=cpostscript -output=$dirname/clara_emit.ps -layout=1,2 -linetypedefault=0,thickness=2 -scales=0,0,-0.2,2.01 -title= -pSpace=0.1,0.95,0.15,0.93 -lspace=0.45,0.52,-0.81,0.49 -graph=line,vary -split=pages -column=s,ecnx -factor=yMultiplier=1e6 -yscale=id=1 $dirname/clara.sig -column=s,Profile -overlay=xmode=norm,yfact=0.04,yoffset=0 -mode=y=linear -graph=line,type=0 $dirname/clara.mag -endpanel -column=s,ecny -scales=0,0,-0.2,4 -factor=yMultiplier=1e6 -yscale=id=1 $dirname/clara.sig -column=s,Profile -overlay=xmode=norm,yfact=0.04,yoffset=0 -mode=y=linear -graph=line,type=0 $dirname/clara.mag >& /dev/null
    
    sddsprocess $dirname/$2.SDDS -process=t,median,medt -process=p,median,medp >& /dev/null
    
    median=`sdds2stream -parameter=medt $dirname/$2.SDDS | tail -n1`
    
    sddsplot -device=cpostscript -output=$dirname/w-init.ps -column=t,p -sparse=1 -graph=dot -filenamesontopline -title= -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=0.511 -xlabel="t (fs)" -ylabel="Energy (MeV)" $dirname/$2.SDDS >& /dev/null
    sddsplot -device=cpostscript -output=$dirname/w-init_xps.ps -column=x,xp -sparse=1 -graph=dot -filenamesontopline -title= $dirname/$2.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/w-init_tx.ps -column=t,x -sparse=1 -graph=dot -filenamesontopline -title= -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e3 -xlabel="t (fs)" -ylabel="x (mm)" $dirname/$2.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/w-init_xy.ps -column=x,y -sparse=1 -graph=dot -filenamesontopline -title= -offset=xChange=-$median -factor=xMultiplier=1e3,yMultiplier=1e3 -xlabel="x (mm)" -ylabel="y (mm)" $dirname/$2.SDDS >& /dev/null
    
    sddshist $dirname/$2.SDDS $dirname/freq_init.SDDS -dataColumn=t -sizeOfBins=20e-15 >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/freq_init.ps -linetypedefault=0,thickness=2 -graph=line -filenamesontopline -column=t,frequency -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1 -xlabel="t (fs)" -ylabel="Current (A)" $dirname/freq_init.SDDS >& /dev/null
    
    sddssort $dirname/$2.SDDS -pipe=out -column=t | sddsprocess -pipe -define=column,bin,"t 20e-15 / int" -process=p,ave,pAve | sddsbreak -pipe -change=bin | sddsprocess -pipe -define=col,delta,"p pAve - pAve /" -process=delta,ave,%sAve | sddsanalyzebeam -pipe=in $dirname/w-init.sliceAnalysis >& /dev/null
    
    sddsprocess $dirname/w-init.sliceAnalysis $dirname/w-init.sliceAnalysis2 -filter=column,Ct,0,1e-8,! >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/slemit_init.ps -linetypedefault=0,thickness=2 -scales=0,0,-0.1,1.5 -lspace=0.52,0.64,-0.05,0.95 -filenamesontopline -graph=line,vary -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e6 -xlabel="t (fs)" -ylabel="\$ge\$r\$bx\,n\$n\, \$ge\$r\$by\,n\$n (mm mrad)" $dirname/w-init.sliceAnalysis2 -column=Ct,enx -legend=specified="\$ge\$r\$bx\,n\$n" -column=Ct,eny -legend=specified="\$ge\$r\$by\,n\$n" >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/slen_init.ps -linetypedefault=0,thickness=2 -filenamesontopline -graph=line,vary -column=Ct,Sdelta -offset=xChange=-$median -factor=xMultiplier=1e15 -xlabel="t (fs)" -ylabel="\$gs\$bd\$n\$r" $dirname/w-init.sliceAnalysis2 >& /dev/null
    
    sddsplot -layout=3,2 -device=cpostscript -output=$dirname/clara_init.ps -pSpace=0.2,0.90,0.2,0.8 -linetypedefault=0,thickness=2 -title= -column=t,p -split=pages -endpanel -scales=0,0,0,0 -sparse=1 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/$2.SDDS -graph=dot,vary -column=t,frequency -split=pages -endpanel -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_init.SDDS -graph=line,vary -column=Ct,\(enx,eny\) -split=pages -endpanel -scales=0,0,-0.1,1.5 -lspace=0.72,0.84,-0.05,0.95 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e6 -xlabel="t (ps)" -ylabel="\$ge\$r\$bx\,n\$n \,\$ge\$r\$by\,n\$n (mm mrad)" $dirname/w-init.sliceAnalysis2 -graph=line,vary -column=Ct,Sdelta -split=pages -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=10000 -xlabel="t (ps)" -ylabel="\$gs\$bd\$n\$r \$sx\$e10\$a4\$n" $dirname/w-init.sliceAnalysis2 -graph=line,vary -endpanel -column=t,x -graph=dot,vary -sparse=1 -split=pages -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e3 $dirname/$2.SDDS -xlabel="t (ps)" -ylabel="x (mm)" -endpanel -column=x,xp -sparse=1 -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e3,yMultiplier=1e3 -xlabel="x (mm)" -ylabel="x' (mrad)" $dirname/$2.SDDS >& /dev/null
	
    sddsprocess $dirname/CLA-S03-DIA-SCR-01.SDDS -process=t,median,medt -process=p,median,medp >& /dev/null
    median=`sdds2stream -parameter=medt $dirname/CLA-S03-DIA-SCR-01.SDDS | tail -n1`

    sddsplot -device=cpostscript -output=$dirname/CLA-S03-DIA-SCR-01.ps -column=t,p -sparse=1 -graph=dot -filenamesontopline -title= -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=0.511 -xlabel="t (fs)" -ylabel="Energy (MeV)" $dirname/CLA-S03-DIA-SCR-01.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-S03-DIA-SCR-01_xps.ps -column=x,xp -sparse=1 -graph=dot -filenamesontopline -title= $dirname/CLA-S03-DIA-SCR-01.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-S03-DIA-SCR-01_tx.ps -column=t,x -sparse=1 -graph=dot -filenamesontopline -title= -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e3 -xlabel="t (fs)" -ylabel="x (mm)" $dirname/CLA-S03-DIA-SCR-01.SDDS >& /dev/null
    
    sddshist $dirname/CLA-S03-DIA-SCR-01.SDDS $dirname/freq_CLA-S03-DIA-SCR-01.SDDS -dataColumn=t -sizeOfBins=20e-15 >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/freq_CLA-S03-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -graph=line -filenamesontopline -column=t,frequency -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1 -xlabel="t (fs)" -ylabel="Current (A)" $dirname/freq_CLA-S03-DIA-SCR-01.SDDS >& /dev/null
    
    sddssort $dirname/CLA-S03-DIA-SCR-01.SDDS -pipe=out -column=t | sddsprocess -pipe -define=column,bin,"t 20e-15 / int" -process=p,ave,pAve | sddsbreak -pipe -change=bin | sddsprocess -pipe -define=col,delta,"p pAve - pAve /" -process=delta,ave,%sAve | sddsanalyzebeam -pipe | sddsbreak -pipe=in $dirname/CLA-S03-DIA-SCR-01.sliceAnalysis -changeOf=Step
    
    sddsprocess $dirname/CLA-S03-DIA-SCR-01.sliceAnalysis $dirname/CLA-S03-DIA-SCR-01.sliceAnalysis2 -filter=column,Ct,0,1e-8,!
    
    sddsplot -device=cpostscript -output=$dirname/slemit_CLA-S03-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -scales=0,0,-0.1,1.5 -lspace=0.52,0.64,-0.05,0.95 -filenamesontopline -graph=line,vary -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e6 -xlabel="t (fs)" -ylabel="\$ge\$r\$bx\,n\$n\, \$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-S03-DIA-SCR-01.sliceAnalysis2 -column=Ct,enx -legend=specified="\$ge\$r\$bx\,n\$n" -column=Ct,eny -legend=specified="\$ge\$r\$by\,n\$n"
    
    sddsplot -device=cpostscript -output=$dirname/slen_CLA-S03-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -filenamesontopline -graph=line,vary -column=Ct,Sdelta -offset=xChange=-$median -factor=xMultiplier=1e15 -xlabel="t (fs)" -ylabel="\$gs\$bd\$n\$r" $dirname/CLA-S03-DIA-SCR-01.sliceAnalysis2
    
    sddsplot -layout=3,2 -device=cpostscript -output=$dirname/clara_CLA-S03-DIA-SCR-01.ps -pSpace=0.2,0.90,0.2,0.8 -linetypedefault=0,thickness=2 -title= -column=t,p -split=pages -endpanel -scales=0,0,0,0 -sparse=1 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-S03-DIA-SCR-01.SDDS -graph=dot,vary -column=t,frequency -split=pages -endpanel -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-S03-DIA-SCR-01.SDDS -graph=line,vary -column=Ct,\(enx,eny\) -split=pages -endpanel -scales=0,0,-0.1,1.5 -lspace=0.72,0.84,-0.05,0.95 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e6 -xlabel="t (ps)" -ylabel="\$ge\$r\$bx\,n\$n \,\$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-S03-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -column=Ct,Sdelta -split=pages -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=10000 -xlabel="t (ps)" -ylabel="\$gs\$bd\$n\$r \$sx\$e10\$a4\$n" $dirname/CLA-S03-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -endpanel -column=dt,x -graph=dot,vary -sparse=1 -split=pages -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e3 $dirname/CLA-S03-DIA-SCR-01.SDDS -xlabel="t (ps)" -ylabel="x (mm)" -endpanel -column=x,xp -sparse=1 -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e3,yMultiplier=1e3 -xlabel="x (mm)" -ylabel="x' (mrad)" $dirname/CLA-S03-DIA-SCR-01.SDDS
	
    sddsprocess $dirname/CLA-S05-DIA-SCR-01.SDDS -process=t,median,medt -process=p,median,medp >& /dev/null
    median=`sdds2stream -parameter=medt $dirname/CLA-S05-DIA-SCR-01.SDDS | tail -n1`

    sddsplot -device=cpostscript -output=$dirname/CLA-S05-DIA-SCR-01.ps -column=t,p -sparse=1 -graph=dot,vary -split=pages -filenamesontopline -title= -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=0.511 -xlabel="t (fs)" -ylabel="Energy (MeV)" $dirname/CLA-S05-DIA-SCR-01.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-S05-DIA-SCR-01_xps.ps -column=x,xp -sparse=1 -graph=dot,vary -split=pages -filenamesontopline -title= $dirname/CLA-S05-DIA-SCR-01.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-S05-DIA-SCR-01_tx.ps -column=t,x -sparse=1 -graph=dot,vary -split=pages -filenamesontopline -title= -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e3 -xlabel="t (fs)" -ylabel="x (mm)" $dirname/CLA-S05-DIA-SCR-01.SDDS >& /dev/null
    
    sddshist $dirname/CLA-S05-DIA-SCR-01.SDDS $dirname/freq_CLA-S05-DIA-SCR-01.SDDS -dataColumn=t -sizeOfBins=20e-15 >& /dev/null
    sddsplot -device=cpostscript -output=$dirname/freq_CLA-S05-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -graph=line,vary -filenamesontopline -column=t,frequency -split=pages -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1 -xlabel="t (fs)" -ylabel="Current (A)" $dirname/freq_CLA-S05-DIA-SCR-01.SDDS >& /dev/null
    
    sddssort $dirname/CLA-S05-DIA-SCR-01.SDDS -pipe=out -column=t | sddsprocess -pipe -define=column,bin,"t 20e-15 / int" -process=p,ave,pAve | sddsbreak -pipe -change=bin | sddsprocess -pipe -define=col,delta,"p pAve - pAve /" -process=delta,ave,%sAve | sddsanalyzebeam -pipe | sddsbreak -pipe=in $dirname/CLA-S05-DIA-SCR-01.sliceAnalysis -changeOf=Step >& /dev/null

    sddsprocess $dirname/CLA-S05-DIA-SCR-01.sliceAnalysis $dirname/CLA-S05-DIA-SCR-01.sliceAnalysis2 -filter=column,Ct,0,1e-8,! >& /dev/null

    sddsplot -device=cpostscript -output=$dirname/slemit_CLA-S05-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -scales=0,0,-0.1,1.5 -lspace=0.52,0.64,-0.05,0.95 -filenamesontopline -graph=line,vary -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e6 -xlabel="t (fs)" -ylabel="\$ge\$r\$bx\,n\$n\, \$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-S05-DIA-SCR-01.sliceAnalysis2 -column=Ct,enx -split=pages -legend=specified="\$ge\$r\$bx\,n\$n" -column=Ct,eny -split=pages -legend=specified="\$ge\$r\$by\,n\$n" >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/slen_CLA-S05-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -filenamesontopline -graph=line,vary -column=Ct,Sdelta -split=pages -offset=xChange=-$median -factor=xMultiplier=1e15 -xlabel="t (fs)" -ylabel="\$gs\$bd\$n\$r" $dirname/CLA-S05-DIA-SCR-01.sliceAnalysis2 >& /dev/null
    
    sddsplot -layout=3,2 -device=cpostscript -output=$dirname/clara_CLA-S05-DIA-SCR-01.ps -pSpace=0.2,0.90,0.2,0.8 -linetypedefault=0,thickness=2 -title= -column=t,p -split=pages -endpanel -scales=0,0,0,0 -sparse=1 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-S05-DIA-SCR-01.SDDS -graph=dot,vary -column=t,frequency -split=pages -endpanel -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-S05-DIA-SCR-01.SDDS -graph=line,vary -column=Ct,\(enx,eny\) -split=pages -endpanel -scales=0,0,-0.1,1.5 -lspace=0.72,0.84,-0.05,0.95 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e6 -xlabel="t (ps)" -ylabel="\$ge\$r\$bx\,n\$n \,\$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-S05-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -column=Ct,Sdelta -split=pages -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=10000 -xlabel="t (ps)" -ylabel="\$gs\$bd\$n\$r \$sx\$e10\$a4\$n" $dirname/CLA-S05-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -endpanel -column=dt,x -graph=dot,vary -sparse=1 -split=pages -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e3 $dirname/CLA-S05-DIA-SCR-01.SDDS -xlabel="t (ps)" -ylabel="x (mm)" -endpanel -column=x,xp -sparse=1 -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e3,yMultiplier=1e3 -xlabel="x (mm)" -ylabel="x' (mrad)" $dirname/CLA-S05-DIA-SCR-01.SDDS >& /dev/null

    sddsprocess $dirname/CLA-S05-DIA-SCR-05.SDDS -process=t,median,medt -process=p,median,medp >& /dev/null
    median=`sdds2stream -parameter=medt $dirname/CLA-S05-DIA-SCR-05.SDDS | tail -n1`

    sddshist $dirname/CLA-S05-DIA-SCR-05.SDDS $dirname/freq_CLA-S05-DIA-SCR-05.SDDS -dataColumn=t -sizeOfBins=20e-15 >& /dev/null
    sddsplot -device=cpostscript -output=$dirname/freq_CLA-S05-DIA-SCR-05.ps -linetypedefault=0,thickness=2 -graph=line,vary -filenamesontopline -column=t,frequency -split=pages -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1 -xlabel="t (fs)" -ylabel="Current (A)" $dirname/freq_CLA-S05-DIA-SCR-05.SDDS >& /dev/null
    
    sddssort $dirname/CLA-S05-DIA-SCR-05.SDDS -pipe=out -column=t | sddsprocess -pipe -define=column,bin,"t 20e-15 / int" -process=p,ave,pAve | sddsbreak -pipe -change=bin | sddsprocess -pipe -define=col,delta,"p pAve - pAve /" -process=delta,ave,%sAve | sddsanalyzebeam -pipe | sddsbreak -pipe=in $dirname/CLA-S05-DIA-SCR-05.sliceAnalysis -changeOf=Step >& /dev/null

    sddsprocess $dirname/CLA-S05-DIA-SCR-05.sliceAnalysis $dirname/CLA-S05-DIA-SCR-05.sliceAnalysis2 -filter=column,Ct,0,1e-8,! >& /dev/null

    sddsplot -layout=3,2 -device=cpostscript -output=$dirname/clara_CLA-S05-DIA-SCR-05.ps -pSpace=0.2,0.90,0.2,0.8 -linetypedefault=0,thickness=2 -title= -column=t,p -split=pages -endpanel -scales=0,0,0,0 -sparse=1 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-S05-DIA-SCR-05.SDDS -graph=dot,vary -column=t,frequency -split=pages -endpanel -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-S05-DIA-SCR-05.SDDS -graph=line,vary -column=Ct,\(enx,eny\) -split=pages -endpanel -scales=0,0,-0.1,1.5 -lspace=0.72,0.84,-0.05,0.95 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e6 -xlabel="t (ps)" -ylabel="\$ge\$r\$bx\,n\$n \,\$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-S05-DIA-SCR-05.sliceAnalysis2 -graph=line,vary -column=Ct,Sdelta -split=pages -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=10000 -xlabel="t (ps)" -ylabel="\$gs\$bd\$n\$r \$sx\$e10\$a4\$n" $dirname/CLA-S05-DIA-SCR-05.sliceAnalysis2 -graph=line,vary -endpanel -column=dt,x -graph=dot,vary -sparse=1 -split=pages -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e3 $dirname/CLA-S05-DIA-SCR-05.SDDS -xlabel="t (ps)" -ylabel="x (mm)" -endpanel -column=x,xp -sparse=1 -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e3,yMultiplier=1e3 -xlabel="x (mm)" -ylabel="x' (mrad)" $dirname/CLA-S05-DIA-SCR-05.SDDS >& /dev/null

    sddsprocess $dirname/CLA-S06-DIA-SCR-01.SDDS -process=t,median,medt -process=p,median,medp >& /dev/null
    median=`sdds2stream -parameter=medt $dirname/CLA-S06-DIA-SCR-01.SDDS | tail -n1`

    sddshist $dirname/CLA-S06-DIA-SCR-01.SDDS $dirname/freq_CLA-S06-DIA-SCR-01.SDDS -dataColumn=t -sizeOfBins=20e-15 >& /dev/null
    sddsplot -device=cpostscript -output=$dirname/freq_CLA-S06-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -graph=line,vary -filenamesontopline -column=t,frequency -split=pages -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1 -xlabel="t (fs)" -ylabel="Current (A)" $dirname/freq_CLA-S06-DIA-SCR-01.SDDS >& /dev/null
    
    sddssort $dirname/CLA-S06-DIA-SCR-01.SDDS -pipe=out -column=t | sddsprocess -pipe -define=column,bin,"t 20e-15 / int" -process=p,ave,pAve | sddsbreak -pipe -change=bin | sddsprocess -pipe -define=col,delta,"p pAve - pAve /" -process=delta,ave,%sAve | sddsanalyzebeam -pipe | sddsbreak -pipe=in $dirname/CLA-S06-DIA-SCR-01.sliceAnalysis -changeOf=Step >& /dev/null

    sddsprocess $dirname/CLA-S06-DIA-SCR-01.sliceAnalysis $dirname/CLA-S06-DIA-SCR-01.sliceAnalysis2 -filter=column,Ct,0,1e-8,! >& /dev/null

    sddsplot -layout=3,2 -device=cpostscript -output=$dirname/clara_CLA-S06-DIA-SCR-01.ps -pSpace=0.2,0.90,0.2,0.8 -linetypedefault=0,thickness=2 -title= -column=t,p -split=pages -endpanel -scales=0,0,0,0 -sparse=1 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-S06-DIA-SCR-01.SDDS -graph=dot,vary -column=t,frequency -split=pages -endpanel -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-S06-DIA-SCR-01.SDDS -graph=line,vary -column=Ct,\(enx,eny\) -split=pages -endpanel -scales=0,0,-0.1,1.5 -lspace=0.72,0.84,-0.05,0.95 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e6 -xlabel="t (ps)" -ylabel="\$ge\$r\$bx\,n\$n \,\$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-S06-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -column=Ct,Sdelta -split=pages -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=10000 -xlabel="t (ps)" -ylabel="\$gs\$bd\$n\$r \$sx\$e10\$a4\$n" $dirname/CLA-S06-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -endpanel -column=dt,x -graph=dot,vary -sparse=1 -split=pages -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e3 $dirname/CLA-S06-DIA-SCR-01.SDDS -xlabel="t (ps)" -ylabel="x (mm)" -endpanel -column=x,xp -sparse=1 -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e3,yMultiplier=1e3 -xlabel="x (mm)" -ylabel="x' (mrad)" $dirname/CLA-S06-DIA-SCR-01.SDDS >& /dev/null

    sddsprocess $dirname/CLA-S07-DIA-SCR-05.SDDS -process=t,median,medt -process=p,median,medp -define=param,phase,"Step 42 / 0.5 - 42 * 5.5 * 11 / 90 +",type=double >& /dev/null
    median=`sdds2stream -parameter=medt $dirname/CLA-S07-DIA-SCR-05.SDDS | tail -n1`
    
    sddsplot -device=cpostscript -output=$dirname/CLA-S07-DIA-SCR-05.ps -column=t,p -split=pages -sparse=1 -graph=dot,vary -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-S07-DIA-SCR-05.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-S07-DIA-SCR-05_xps.ps -column=x,xp -sparse=1 -graph=dot,vary -split=pages -filenamesontopline -title= $dirname/CLA-S07-DIA-SCR-05.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-S07-DIA-SCR-05_tx.ps -column=t,x -split=pages -sparse=1 -graph=dot,vary -filenamesontopline -title= -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e3 -xlabel="t (fs)" -ylabel="x (mm)" $dirname/CLA-S07-DIA-SCR-05.SDDS >& /dev/null
    
    sddshist $dirname/CLA-S07-DIA-SCR-05.SDDS $dirname/freq_CLA-S07-DIA-SCR-05.SDDS -dataColumn=t -sizeOfBins=20e-15 >& /dev/null
    sddsplot -device=cpostscript -output=$dirname/freq_CLA-S07-DIA-SCR-05.ps -linetypedefault=0,thickness=2 -graph=line,vary -split=pages -filenamesontopline -column=t,frequency -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-S07-DIA-SCR-05.SDDS >& /dev/null
    
    sddssort $dirname/CLA-S07-DIA-SCR-05.SDDS -pipe=out -column=t | sddsprocess -pipe -define=column,bin,"t 20e-15 / int" -process=p,ave,pAve | sddsbreak -pipe -change=bin | sddsprocess -pipe -define=col,delta,"p pAve - pAve /" -process=delta,ave,%sAve | sddsanalyzebeam -pipe | sddsbreak -pipe=in $dirname/CLA-S07-DIA-SCR-05.sliceAnalysis -changeOf=Step >& /dev/null

    sddsprocess $dirname/CLA-S07-DIA-SCR-05.sliceAnalysis $dirname/CLA-S07-DIA-SCR-05.sliceAnalysis2 -filter=column,Ct,0,1e-8,! >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/slemit_CLA-S07-DIA-SCR-05.ps -linetypedefault=0,thickness=2 -scales=0,0,-0.1,1.5 -lspace=0.52,0.64,-0.05,0.95 -filenamesontopline -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e6 -xlabel="t (fs)" -ylabel="\$ge\$r\$bx\,n\$n\, \$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-S07-DIA-SCR-05.sliceAnalysis2 -column=Ct,enx -graph=line,vary -split=pages -legend=specified="\$ge\$r\$bx\,n\$n" -column=Ct,eny -graph=line,vary -split=pages -legend=specified="\$ge\$r\$by\,n\$n" >& /dev/null

    sddsplot -device=cpostscript -output=$dirname/slen_CLA-S07-DIA-SCR-05.ps -linetypedefault=0,thickness=2 -filenamesontopline -column=Ct,Sdelta -graph=line,vary -split=pages -offset=xChange=-$median -factor=xMultiplier=1e15 -xlabel="t (fs)" -ylabel="\$gs\$bd\$n\$r" $dirname/CLA-S07-DIA-SCR-05.sliceAnalysis2 >& /dev/null
    
    sddsplot -layout=3,2 -device=cpostscript -output=$dirname/clara_CLA-S07-DIA-SCR-05.ps -pSpace=0.2,0.90,0.2,0.8 -linetypedefault=0,thickness=2 -title= -column=t,p -split=pages -endpanel -scales=0,0,0,0 -sparse=1 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-S07-DIA-SCR-05.SDDS -graph=dot,vary -column=t,frequency -split=pages -endpanel -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-S07-DIA-SCR-05.SDDS -graph=line,vary -column=Ct,\(enx,eny\) -split=pages -endpanel -scales=0,0,-0.1,1.5 -lspace=0.72,0.84,-0.05,0.95 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e6 -xlabel="t (ps)" -ylabel="\$ge\$r\$bx\,n\$n \,\$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-S07-DIA-SCR-05.sliceAnalysis2 -graph=line,vary -column=Ct,Sdelta -split=pages -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=10000 -xlabel="t (ps)" -ylabel="\$gs\$bd\$n\$r \$sx\$e10\$a4\$n" $dirname/CLA-S07-DIA-SCR-05.sliceAnalysis2 -graph=line,vary -endpanel -column=dt,x -graph=dot,vary -sparse=1 -split=pages -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e3 $dirname/CLA-S07-DIA-SCR-05.SDDS -xlabel="t (ps)" -ylabel="x (mm)" -endpanel -column=x,xp -sparse=1 -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e3,yMultiplier=1e3 -xlabel="x (mm)" -ylabel="x' (mrad)" $dirname/CLA-S07-DIA-SCR-05.SDDS >& /dev/null

    sddsprocess $dirname/CLA-MU1-DIA-SCR-01.SDDS -process=t,median,medt -process=p,median,medp -define=param,phase,"Step 42 / 0.5 - 42 * 5.5 * 11 / 90 +",type=double >& /dev/null
    median=`sdds2stream -parameter=medt $dirname/CLA-MU1-DIA-SCR-01.SDDS | tail -n1`
    
    sddsplot -device=cpostscript -output=$dirname/CLA-MU1-DIA-SCR-01.ps -column=t,p -split=pages -sparse=1 -graph=dot,vary -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-MU1-DIA-SCR-01.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-MU1-DIA-SCR-01_xps.ps -column=x,xp -sparse=1 -graph=dot,vary -split=pages -filenamesontopline -title= $dirname/CLA-MU1-DIA-SCR-01.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-MU1-DIA-SCR-01_tx.ps -column=t,x -split=pages -sparse=1 -graph=dot,vary -filenamesontopline -title= -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e3 -xlabel="t (fs)" -ylabel="x (mm)" $dirname/CLA-MU1-DIA-SCR-01.SDDS >& /dev/null
    
    sddshist $dirname/CLA-MU1-DIA-SCR-01.SDDS $dirname/freq_CLA-MU1-DIA-SCR-01.SDDS -dataColumn=t -sizeOfBins=20e-15 >& /dev/null
    sddsplot -device=cpostscript -output=$dirname/freq_CLA-MU1-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -graph=line,vary -split=pages -filenamesontopline -column=t,frequency -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-MU1-DIA-SCR-01.SDDS >& /dev/null
    
    sddssort $dirname/CLA-MU1-DIA-SCR-01.SDDS -pipe=out -column=t | sddsprocess -pipe -define=column,bin,"t 20e-15 / int" -process=p,ave,pAve | sddsbreak -pipe -change=bin | sddsprocess -pipe -define=col,delta,"p pAve - pAve /" -process=delta,ave,%sAve | sddsanalyzebeam -pipe | sddsbreak -pipe=in $dirname/CLA-MU1-DIA-SCR-01.sliceAnalysis -changeOf=Step >& /dev/null

    sddsprocess $dirname/CLA-MU1-DIA-SCR-01.sliceAnalysis $dirname/CLA-MU1-DIA-SCR-01.sliceAnalysis2 -filter=column,Ct,0,1e-8,! >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/slemit_CLA-MU1-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -scales=0,0,-0.1,1.5 -lspace=0.52,0.64,-0.05,0.95 -filenamesontopline -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e6 -xlabel="t (fs)" -ylabel="\$ge\$r\$bx\,n\$n\, \$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-MU1-DIA-SCR-01.sliceAnalysis2 -column=Ct,enx -graph=line,vary -split=pages -legend=specified="\$ge\$r\$bx\,n\$n" -column=Ct,eny -graph=line,vary -split=pages -legend=specified="\$ge\$r\$by\,n\$n" >& /dev/null

    sddsplot -device=cpostscript -output=$dirname/slen_CLA-MU1-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -filenamesontopline -column=Ct,Sdelta -graph=line,vary -split=pages -offset=xChange=-$median -factor=xMultiplier=1e15 -xlabel="t (fs)" -ylabel="\$gs\$bd\$n\$r" $dirname/CLA-MU1-DIA-SCR-01.sliceAnalysis2 >& /dev/null
    
    sddsplot -layout=3,2 -device=cpostscript -output=$dirname/clara_CLA-MU1-DIA-SCR-01.ps -pSpace=0.2,0.90,0.2,0.8 -linetypedefault=0,thickness=2 -title= -column=t,p -split=pages -endpanel -scales=0,0,0,0 -sparse=1 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-MU1-DIA-SCR-01.SDDS -graph=dot,vary -column=t,frequency -split=pages -endpanel -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-MU1-DIA-SCR-01.SDDS -graph=line,vary -column=Ct,\(enx,eny\) -split=pages -endpanel -scales=0,0,-0.1,1.5 -lspace=0.72,0.84,-0.05,0.95 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e6 -xlabel="t (ps)" -ylabel="\$ge\$r\$bx\,n\$n \,\$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-MU1-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -column=Ct,Sdelta -split=pages -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=10000 -xlabel="t (ps)" -ylabel="\$gs\$bd\$n\$r \$sx\$e10\$a4\$n" $dirname/CLA-MU1-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -endpanel -column=dt,x -graph=dot,vary -sparse=1 -split=pages -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e3 $dirname/CLA-MU1-DIA-SCR-01.SDDS -xlabel="t (ps)" -ylabel="x (mm)" -endpanel -column=x,xp -sparse=1 -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e3,yMultiplier=1e3 -xlabel="x (mm)" -ylabel="x' (mrad)" $dirname/CLA-MU1-DIA-SCR-01.SDDS >& /dev/null

    sddsprocess $dirname/CLA-RU1-DIA-SCR-01.SDDS -process=t,median,medt -process=p,median,medp -define=param,phase,"Step 42 / 0.5 - 42 * 5.5 * 11 / 90 +",type=double >& /dev/null
    median=`sdds2stream -parameter=medt $dirname/CLA-RU1-DIA-SCR-01.SDDS | tail -n1`
    
    sddsplot -device=cpostscript -output=$dirname/CLA-RU1-DIA-SCR-01.ps -column=t,p -split=pages -sparse=1 -graph=dot,vary -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-RU1-DIA-SCR-01.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-RU1-DIA-SCR-01_xps.ps -column=x,xp -sparse=1 -graph=dot,vary -split=pages -filenamesontopline -title= $dirname/CLA-RU1-DIA-SCR-01.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-RU1-DIA-SCR-01_tx.ps -column=t,x -split=pages -sparse=1 -graph=dot,vary -filenamesontopline -title= -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e3 -xlabel="t (fs)" -ylabel="x (mm)" $dirname/CLA-RU1-DIA-SCR-01.SDDS >& /dev/null
    
    sddshist $dirname/CLA-RU1-DIA-SCR-01.SDDS $dirname/freq_CLA-RU1-DIA-SCR-01.SDDS -dataColumn=t -sizeOfBins=20e-15 >& /dev/null
    sddsplot -device=cpostscript -output=$dirname/freq_CLA-RU1-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -graph=line,vary -split=pages -filenamesontopline -column=t,frequency -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-RU1-DIA-SCR-01.SDDS >& /dev/null
    
    sddssort $dirname/CLA-RU1-DIA-SCR-01.SDDS -pipe=out -column=t | sddsprocess -pipe -define=column,bin,"t 20e-15 / int" -process=p,ave,pAve | sddsbreak -pipe -change=bin | sddsprocess -pipe -define=col,delta,"p pAve - pAve /" -process=delta,ave,%sAve | sddsanalyzebeam -pipe | sddsbreak -pipe=in $dirname/CLA-RU1-DIA-SCR-01.sliceAnalysis -changeOf=Step >& /dev/null

    sddsprocess $dirname/CLA-RU1-DIA-SCR-01.sliceAnalysis $dirname/CLA-RU1-DIA-SCR-01.sliceAnalysis2 -filter=column,Ct,0,1e-8,! >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/slemit_CLA-RU1-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -scales=0,0,-0.1,1.5 -lspace=0.52,0.64,-0.05,0.95 -filenamesontopline -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e6 -xlabel="t (fs)" -ylabel="\$ge\$r\$bx\,n\$n\, \$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-RU1-DIA-SCR-01.sliceAnalysis2 -column=Ct,enx -graph=line,vary -split=pages -legend=specified="\$ge\$r\$bx\,n\$n" -column=Ct,eny -graph=line,vary -split=pages -legend=specified="\$ge\$r\$by\,n\$n" >& /dev/null

    sddsplot -device=cpostscript -output=$dirname/slen_CLA-RU1-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -filenamesontopline -column=Ct,Sdelta -graph=line,vary -split=pages -offset=xChange=-$median -factor=xMultiplier=1e15 -xlabel="t (fs)" -ylabel="\$gs\$bd\$n\$r" $dirname/CLA-RU1-DIA-SCR-01.sliceAnalysis2 >& /dev/null
    
    sddsplot -layout=3,2 -device=cpostscript -output=$dirname/clara_CLA-RU1-DIA-SCR-01.ps -pSpace=0.2,0.90,0.2,0.8 -linetypedefault=0,thickness=2 -title= -column=t,p -split=pages -endpanel -scales=0,0,0,0 -sparse=1 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-RU1-DIA-SCR-01.SDDS -graph=dot,vary -column=t,frequency -split=pages -endpanel -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-RU1-DIA-SCR-01.SDDS -graph=line,vary -column=Ct,\(enx,eny\) -split=pages -endpanel -scales=0,0,-0.1,1.5 -lspace=0.72,0.84,-0.05,0.95 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e6 -xlabel="t (ps)" -ylabel="\$ge\$r\$bx\,n\$n \,\$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-RU1-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -column=Ct,Sdelta -split=pages -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=10000 -xlabel="t (ps)" -ylabel="\$gs\$bd\$n\$r \$sx\$e10\$a4\$n" $dirname/CLA-RU1-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -endpanel -column=dt,x -graph=dot,vary -sparse=1 -split=pages -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e3 $dirname/CLA-RU1-DIA-SCR-01.SDDS -xlabel="t (ps)" -ylabel="x (mm)" -endpanel -column=x,xp -sparse=1 -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e3,yMultiplier=1e3 -xlabel="x (mm)" -ylabel="x' (mrad)" $dirname/CLA-RU1-DIA-SCR-01.SDDS >& /dev/null

    sddsprocess $dirname/CLA-RU7-DIA-SCR-01.SDDS -process=t,median,medt -process=p,median,medp -define=param,phase,"Step 42 / 0.5 - 42 * 5.5 * 11 / 90 +",type=double >& /dev/null
    median=`sdds2stream -parameter=medt $dirname/CLA-RU7-DIA-SCR-01.SDDS | tail -n1`
    
    sddsplot -device=cpostscript -output=$dirname/CLA-RU7-DIA-SCR-01.ps -column=t,p -split=pages -sparse=1 -graph=dot,vary -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-RU7-DIA-SCR-01.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-RU7-DIA-SCR-01_xps.ps -column=x,xp -sparse=1 -graph=dot,vary -split=pages -filenamesontopline -title= $dirname/CLA-RU7-DIA-SCR-01.SDDS >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/CLA-RU7-DIA-SCR-01_tx.ps -column=t,x -split=pages -sparse=1 -graph=dot,vary -filenamesontopline -title= -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e3 -xlabel="t (fs)" -ylabel="x (mm)" $dirname/CLA-RU7-DIA-SCR-01.SDDS >& /dev/null
    
    sddshist $dirname/CLA-RU7-DIA-SCR-01.SDDS $dirname/freq_CLA-RU7-DIA-SCR-01.SDDS -dataColumn=t -sizeOfBins=20e-15 >& /dev/null
    sddsplot -device=cpostscript -output=$dirname/freq_CLA-RU7-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -graph=line,vary -split=pages -filenamesontopline -column=t,frequency -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-RU7-DIA-SCR-01.SDDS >& /dev/null
    
    sddssort $dirname/CLA-RU7-DIA-SCR-01.SDDS -pipe=out -column=t | sddsprocess -pipe -define=column,bin,"t 20e-15 / int" -process=p,ave,pAve | sddsbreak -pipe -change=bin | sddsprocess -pipe -define=col,delta,"p pAve - pAve /" -process=delta,ave,%sAve | sddsanalyzebeam -pipe | sddsbreak -pipe=in $dirname/CLA-RU7-DIA-SCR-01.sliceAnalysis -changeOf=Step >& /dev/null

    sddsprocess $dirname/CLA-RU7-DIA-SCR-01.sliceAnalysis $dirname/CLA-RU7-DIA-SCR-01.sliceAnalysis2 -filter=column,Ct,0,1e-8,! >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/slemit_CLA-RU7-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -scales=0,0,-0.1,1.5 -lspace=0.52,0.64,-0.05,0.95 -filenamesontopline -offset=xChange=-$median -factor=xMultiplier=1e15,yMultiplier=1e6 -xlabel="t (fs)" -ylabel="\$ge\$r\$bx\,n\$n\, \$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-RU7-DIA-SCR-01.sliceAnalysis2 -column=Ct,enx -graph=line,vary -split=pages -legend=specified="\$ge\$r\$bx\,n\$n" -column=Ct,eny -graph=line,vary -split=pages -legend=specified="\$ge\$r\$by\,n\$n" >& /dev/null

    sddsplot -device=cpostscript -output=$dirname/slen_CLA-RU7-DIA-SCR-01.ps -linetypedefault=0,thickness=2 -filenamesontopline -column=Ct,Sdelta -graph=line,vary -split=pages -offset=xChange=-$median -factor=xMultiplier=1e15 -xlabel="t (fs)" -ylabel="\$gs\$bd\$n\$r" $dirname/CLA-RU7-DIA-SCR-01.sliceAnalysis2 >& /dev/null
    
    sddsplot -layout=3,2 -device=cpostscript -output=$dirname/clara_CLA-RU7-DIA-SCR-01.ps -pSpace=0.2,0.90,0.2,0.8 -linetypedefault=0,thickness=2 -title= -column=t,p -split=pages -endpanel -scales=0,0,0,0 -sparse=1 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=0.511 -xlabel="t (ps)" -ylabel="Energy (MeV)" $dirname/CLA-RU7-DIA-SCR-01.SDDS -graph=dot,vary -column=t,frequency -split=pages -endpanel -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1 -xlabel="t (ps)" -ylabel="Current (A)" $dirname/freq_CLA-RU7-DIA-SCR-01.SDDS -graph=line,vary -column=Ct,\(enx,eny\) -split=pages -endpanel -scales=0,0,-0.1,1.5 -lspace=0.72,0.84,-0.05,0.95 -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e6 -xlabel="t (ps)" -ylabel="\$ge\$r\$bx\,n\$n \,\$ge\$r\$by\,n\$n (mm mrad)" $dirname/CLA-RU7-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -column=Ct,Sdelta -split=pages -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=10000 -xlabel="t (ps)" -ylabel="\$gs\$bd\$n\$r \$sx\$e10\$a4\$n" $dirname/CLA-RU7-DIA-SCR-01.sliceAnalysis2 -graph=line,vary -endpanel -column=dt,x -graph=dot,vary -sparse=1 -split=pages -title= -offset=xChange=-$median -factor=xMultiplier=1e12,yMultiplier=1e3 $dirname/CLA-RU7-DIA-SCR-01.SDDS -xlabel="t (ps)" -ylabel="x (mm)" -endpanel -column=x,xp -sparse=1 -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e3,yMultiplier=1e3 -xlabel="x (mm)" -ylabel="x' (mrad)" $dirname/CLA-RU7-DIA-SCR-01.SDDS >& /dev/null
	
    sddsprocess $dirname/$2.SDDS $dirname/w-init.ddt.SDDS -redefine=column,dt,"t medt -" -define=column,dp,"p medp - medp /" >& /dev/null
    sddsprocess $dirname/CLA-S03-DIA-SCR-01.SDDS $dirname/CLA-S03-DIA-SCR-01.ddt.SDDS -redefine=column,dt,"dt medt -" -define=column,dp,"p medp - medp /" >& /dev/null
    sddsprocess $dirname/CLA-S05-DIA-SCR-01.SDDS $dirname/CLA-S05-DIA-SCR-01.ddt.SDDS -redefine=column,dt,"dt medt -" -define=column,dp,"p medp - medp /" >& /dev/null
    sddsprocess $dirname/CLA-S07-DIA-SCR-05.SDDS $dirname/CLA-S07-DIA-SCR-05.ddt.SDDS -redefine=column,dt,"dt medt -" -define=column,dp,"p medp - medp /" >& /dev/null
    sddsprocess $dirname/CLA-MU1-DIA-SCR-01.SDDS $dirname/CLA-MU1-DIA-SCR-01.ddt.SDDS -redefine=column,dt,"dt medt -" -define=column,dp,"p medp - medp /" >& /dev/null
    sddsprocess $dirname/CLA-RU1-DIA-SCR-01.SDDS $dirname/CLA-RU1-DIA-SCR-01.ddt.SDDS -redefine=column,dt,"dt medt -" -define=column,dp,"p medp - medp /" >& /dev/null
    sddsprocess $dirname/CLA-RU7-DIA-SCR-01.SDDS $dirname/CLA-RU7-DIA-SCR-01.ddt.SDDS -redefine=column,dt,"dt medt -" -define=column,dp,"p medp - medp /" >& /dev/null

    sddsplot -device=cpostscript -output=$dirname/clara_lps.ps -sparse=1 -column=dt,p -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e12,yMultiplier=0.511 $dirname/w-init.ddt.SDDS $dirname/CLA-S03-DIA-SCR-01.ddt.SDDS $dirname/CLA-S05-DIA-SCR-01.ddt.SDDS $dirname/CLA-S07-DIA-SCR-05.ddt.SDDS $dirname/CLA-MU1-DIA-SCR-01.ddt.SDDS $dirname/CLA-RU1-DIA-SCR-01.ddt.SDDS $dirname/CLA-RU7-DIA-SCR-01.ddt.SDDS -xlabel="t (ps)" -ylabel="Energy (MeV)" >& /dev/null

    sddsplot -device=cpostscript -output=$dirname/clara_lps_dp.ps -sparse=1 -column=dt,dp -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e12,yMultiplier=100 $dirname/w-init.ddt.SDDS $dirname/CLA-S03-DIA-SCR-01.ddt.SDDS $dirname/CLA-S05-DIA-SCR-01.ddt.SDDS $dirname/CLA-S07-DIA-SCR-05.ddt.SDDS $dirname/CLA-MU1-DIA-SCR-01.ddt.SDDS $dirname/CLA-RU1-DIA-SCR-01.ddt.SDDS $dirname/CLA-RU7-DIA-SCR-01.ddt.SDDS -xlabel="t (ps)" -ylabel="\$gd\$rp (%)" >& /dev/null

    sddsplot -device=cpostscript -output=$dirname/clara_end_dp.ps -sparse=1 -column=dt,dp -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e12,yMultiplier=100 $dirname/CLA-S07-DIA-SCR-05.ddt.SDDS -xlabel="t (ps)" -ylabel="\$gd\$rp (%)" >& /dev/null

    sddsplot -device=cpostscript -output=$dirname/clara_tx.ps -sparse=1 -column=dt,x -graph=dot,vary -split=pages -title= -factor=xMultiplier=1e12,yMultiplier=1e3 $dirname/w-init.ddt.SDDS $dirname/CLA-S03-DIA-SCR-01.ddt.SDDS $dirname/CLA-S05-DIA-SCR-01.ddt.SDDS $dirname/CLA-S07-DIA-SCR-05.ddt.SDDS $dirname/CLA-MU1-DIA-SCR-01.ddt.SDDS $dirname/CLA-RU1-DIA-SCR-01.ddt.SDDS $dirname/CLA-RU7-DIA-SCR-01.ddt.SDDS -xlabel="t (ps)" -ylabel="x (mm)" >& /dev/null
    
    sddsplot -layout=2,2 -device=cpostscript -output=$dirname/clara_bucp_dip.ps -linetypedefault=0,thickness=1.5 -toptitle -col=s,LinearDensity -endpanel -graph=line,vary -split=pages -order=spectral -title="Dipole 1" $dirname/CLA-VBC-MAG-DIP-01.CSR -factor=xMultiplier=1e3 -xlabel="s (mm)" -ylabel="Current (A)" -col=s,LinearDensity -endpanel -graph=line,vary -split=pages -order=spectral -title="Dipole 2" $dirname/CLA-VBC-MAG-DIP-02.CSR -factor=xMultiplier=1e3 -xlabel="s (mm)" -ylabel="Current (A)" -col=s,LinearDensity -endpanel -graph=line,vary -split=pages -order=spectral -title="Dipole 3" $dirname/CLA-VBC-MAG-DIP-03.CSR -factor=xMultiplier=1e3 -xlabel="s (mm)" -ylabel="Current (A)" -col=s,LinearDensity -graph=line,vary -split=pages -order=spectral -title="Dipole 4" $dirname/CLA-VBC-MAG-DIP-04.CSR -factor=xMultiplier=1e3 -xlabel="s (mm)" -ylabel="Current (A)" >& /dev/null
    
    sddsplot -layout=2,2 -device=cpostscript -output=$dirname/clara_bucp_dip2.ps -linetypedefault=0,thickness=1.5 -toptitle -col=s,DeltaGamma -endpanel -graph=line,vary -split=pages -order=spectral -title="Dipole 1" $dirname/CLA-VBC-MAG-DIP-01.CSR -factor=xMultiplier=1e3 -xlabel="s (mm)" -ylabel="\$gDg\$r" -col=s,DeltaGamma -endpanel -graph=line,vary -split=pages -order=spectral -title="Dipole 2" $dirname/CLA-VBC-MAG-DIP-02.CSR -factor=xMultiplier=1e3 -xlabel="s (mm)" -ylabel="\$gDg\$r" -col=s,DeltaGamma -endpanel -graph=line,vary -split=pages -order=spectral -title="Dipole 3" $dirname/CLA-VBC-MAG-DIP-03.CSR -factor=xMultiplier=1e3 -xlabel="s (mm)" -ylabel="\$gDg\$r" -col=s,DeltaGamma -graph=line,vary -split=pages -order=spectral -title="Dipole 4" $dirname/CLA-VBC-MAG-DIP-04.CSR -factor=xMultiplier=1e3 -xlabel="s (mm)" -ylabel="\$gDg\$r"  >& /dev/null
    
    rm -f $dirname/clara_l1.dat 
    touch $dirname/clara_l1.dat
    for ((i=1;i<=9;i++))
      do
      s1=$i
      s2=`sdds2stream -param=gainLengthSlice0$i $dirname/clara.xie`
      echo $s1 $s2 >> $dirname/clara_l1.dat
    done
    for ((i=10;i<=99;i++))
      do
      s1=$i
      s2=`sdds2stream -param=gainLengthSlice$i $dirname/clara.xie`
      echo $s1 $s2 >> $dirname/clara_l1.dat
    done
    plaindata2sdds $dirname/clara_l1.dat $dirname/clara_l1.SDDS -norowcount -column=slice,long -column=gainLength,double
    
    sddsplot -layout=3,2 -device=cpostscript -output=$dirname/clara_xie.ps -pSpace=0.2,0.90,0.2,0.8 -linetypedefault=0,thickness=2 -title= -col=slice,gainLength $dirname/clara_l1.SDDS -ylabel="Gain Length (m)" >& /dev/null
    
    sddsxref $dirname/clara.slan $dirname/clara.cen $dirname/clara.slan2 -take=pCentral -match=ElementName -nowarning >& /dev/null
    
    sddsprocess $dirname/clara.slan2 $dirname/clara.slan3 -define=column,currentSlice%s,"chargeSlice%s durationSlice%s /",select=chargeSlice*,edit=%/chargeSlice//,symbol="Current in Slice %s (A)" -process=enx,first,initemit >& /dev/null
    
    sddsprocess $dirname/clara.slan3 $dirname/clara.slan4 -define=column,laminarity%s,"currentSlice%s 2 / 17000 / initemit / pCentral / 30 / 0.5 / sqr",select=currentSlice*,edit=%/currentSlice//,symbol="Laminarity of Slice %s" >& /dev/null
    
    rm -f $dirname/clara_lam.dat 
    touch $dirname/clara_lam.dat
    for ((i=1;i<=9;i++))
      do
      s1="0"$i
      s2=`sdds2stream -col=currentSlice$s1 $dirname/clara.slan4 | tail -n1`
      echo $s1 $s2 >> $dirname/clara_lam.dat
    done
    i=10
    s1=$i
    s2=`sdds2stream -col=currentSlice$i $dirname/clara.slan4 | tail -n1`
    echo $s1 $s2 >> $dirname/clara_lam.dat
    plaindata2sdds $dirname/clara_lam.dat $dirname/clara_lam.SDDS -norowcount -column=slice,string -column=currentSlice,double
    sddsprocess $dirname/clara_lam.SDDS -process=currentSlice,maximum,maxcurrentslice,functionOf=slice,position -nowarnings
    peakslice=`sdds2stream -param=maxcurrentslice $dirname/clara_lam.SDDS`
    
    sddsplot -device=cpostscript -output=$dirname/clara_lam_peak.ps -lspace=0.6,0.92,-0.91,0.79 -layout=1,2 -linetypedefault=0,thickness=2 -title= -pSpace=0.15,0.93,0.15,0.93 -graph=line,vary  -column=s,laminarity$peakslice -ylabel="\$gr\$r\$bL\$n" $dirname/clara.slan4 -column=s,currentSlice$peakslice -yScalesGroup=id=current $dirname/clara.slan4 -unsuppresszero -ylabel="I\$bPeak\$n (A)" -column=s,Profile -overlay=xmode=norm,yfact=0.04,yoffset=0 -graph=line,type=0 $dirname/clara.mag -scales=0,0,0,0 >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/clara_lam_all.ps -lspace=0.1,0.32,-0.36,0.94 -layout=1,2 -linetypedefault=0,thickness=2 -title= -pSpace=0.15,0.93,0.15,0.93 -graph=line,vary  -column=s,laminarity* -endpanel -legend -ylabel="\$gr\$r\$bL\$n" $dirname/clara.slan4 -column=s,currentSlice* $dirname/clara.slan4 -unsuppresszero -ylabel="I (A)" -legend -column=s,Profile -overlay=xmode=norm,yfact=0.04,yoffset=0 -graph=line,type=0 $dirname/clara.mag -drawline=x0Value=0,x1Value=23,y0Value=1,y1Value=1,linetype=2,thickness=1 -scales=0,0,-40,1000 >& /dev/null
    
    # This makes the gradient vs energy plot - only works for cav scan!
    sddsprocess $dirname/clara.params -pipe=out -define=column,Step,i_page,type=long | sddsprocess -pipe -match=column,ElementName="LIN2-CAV" -match=column,ElementParameter="VOLT" | sddsregroup -pipe=in $dirname/clara.params2  >& /dev/null
    sddscollapse $dirname/CLA-S07-DIA-SCR-05.SDDS $dirname/CLA-S07-DIA-SCR-05_2.SDDS  >& /dev/null
    sddsxref $dirname/CLA-S07-DIA-SCR-05_2.SDDS $dirname/clara.params2 $dirname/clara.escan -take=ParameterValue >& /dev/null
    sddsprocess $dirname/clara.escan -define=column,Gradient,"ParameterValue 0.0333333 / 122 /"  >& /dev/null

    sddsplot -device=cpostscript -output=$dirname/clara_gradient.ps -lspace=0.6,0.92,-0.91,0.79 -layout=1,1 -linetypedefault=0,thickness=2 -column=Gradient,medp -graph=line -title= -ticksettings=grid -factor=xMultiplier=1e-6,yMultiplier=0.511 $dirname/clara.escan -xlabel="Linac 3/4 Gradient (MV/m)" -ylabel="Energy (MeV)" >& /dev/null
    
    sddsplot -device=cpostscript -output=$dirname/clara_lost.ps -layout=1,2 -linetypedefault=0,thickness=2 -scales=0,0,-0.2,2.01 -title= -pSpace=0.1,0.95,0.15,0.93 -lspace=0.45,0.52,-0.81,0.49 -graph=line,vary -split=pages -column=s,x -factor=yMultiplier=1e3 -yscale=id=1 -ylabel="x (mm)" $dirname/clara.lost -column=s,Profile -overlay=xmode=norm,yfact=0.04,yoffset=0 -mode=y=linear -graph=line,type=0 $dirname/clara.mag -endpanel -column=s,y -scales=0,0,-0.2,4 -factor=yMultiplier=1e3 -yscale=id=1 -ylabel="y (mm)" $dirname/clara.lost -column=s,Profile -overlay=xmode=norm,yfact=0.04,yoffset=0 -mode=y=linear -graph=line,type=0 $dirname/clara.mag >& /dev/null


    echo "Analysis in directory $dirname finished"
else
    echo "Elegant did not finish"
fi
if [ $png -eq 1 ];then
    echo "Converting .ps to .png"
    for INP in $dirname/*.ps
      do
      newname=`basename $INP .ps`
      convert -rotate '90<' -density 200 -geometry 100% $INP $dirname/$newname.png
      echo "Converted $INP to $dirname/$newname.png"
    done
    echo "Converted all .ps to .png"
fi
