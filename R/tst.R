


#fa30 <- Sys.glob("/wrk/yuanzhen/figureone/*30n.gz")
#fa101 <- Sys.glob(paste0("/wrk/yuanzhen/figureone/*101n.gz"))

#fjComm::clear_()
#for(p in 1:length(fa30)){
  #seqs=fjComm::getSeq_fqfachrFile(fa30[p])
  #seqs=length_adjust(seqs,output_length = 30,seq_in_middle = TRUE)
  #E_MI=fjComm::ic_related_calc(seqs = seqs,kmerLen = 3L,filter_for_spacing =F ,pseudo=10, type = "maxBias")
  #EMI_plot=gg_heat2D_MI_rotate(E_MI,grad_colors =c("white","steelblue","red"))+guides(fill="none")
  #ggsave(filename=paste0(p,"_30n",".png"),plot=EMI_plot,path="/wrk/yuanzhen/testfg/30pic/")
#}

#file="/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned/B8_20_Mock__root__conc_0nM__time_2h__101n.gz"
#seqs=fjComm::getSeq_fqfachrFile(fa30)

#length_adjust<-function(seqs,output_length=40L,seq_in_middle=FALSE)
  # adjust all ligands to the same length, out put seq only
#{
  #seqs %<>% str_sub(1,output_length)
  #Ns=output_length-nchar(seqs)
  #topaste=strrep("N",Ns)
  #if(seq_in_middle){
    #leftNs=(Ns/2) %>% as.integer()
   # rightNs=Ns-leftNs
    #leftNs=strrep("N",leftNs)
    #rightNs=strrep("N",rightNs)
    #seqs=paste0(leftNs,seqs,rightNs)
  #}else seqs=paste0(seqs,topaste)
  #seqs
#}

#seqs=length_adjust(seqs,output_length = 101,seq_in_middle = TRUE)


#test=F

#color_mid=c("#00700b","steelblue","#838100") %>% set_names(c("MI","maxBias","dimer"))

#total_MI=fjComm::ic_related_calc(seqs = seqs,kmerLen = 3L,filter_for_spacing =ifelse(test,T,F) ,spacing = 0:10,
                                 #verbose=F, pseudo=10, type = "MI")
#MI_plot=gg_heat2D_MI_rotate(total_MI,grad_colors =c("white",color_mid["MI"],"red"))+guides(fill="none")
# gg_heat2D_MI(total_MI)

#E_MI=fjComm::ic_related_calc(seqs = seqs,kmerLen = 3L,filter_for_spacing =ifelse(test,T,F) ,spacing = 0:10,
                             #verbose=F, pseudo=10, type = "maxBias", maxBias_dimer_Params=list(type="topMI",topNo=10L))
#EMI_plot=gg_heat2D_MI_rotate(E_MI,grad_colors =c("white",color_mid["maxBias"],"red"))+guides(fill="none")
#
#dimer_MI=fjComm::ic_related_calc(seqs = seqs,kmerLen = 3L,filter_for_spacing =ifelse(test,T,F) ,spacing = 0:10,
                                 #verbose=F, pseudo=10, type = "dimer", maxBias_dimer_Params=list(type="topMI",topNo=10L))
#dimer_plot=gg_heat2D_MI_rotate(total_MI,grad_colors =c("white",color_mid["dimer"],"red"))+guides(fill="none")

# ggsave(filename = "dimerMIplot.pdf",plot = dimer_plot)


#simplified version
#E_MI=fjComm::ic_related_calc(seqs = seqs,kmerLen = 3L,filter_for_spacing =F ,pseudo=10, type = "maxBias")
#EMI_plot=gg_heat2D_MI_rotate(E_MI,grad_colors =c("white","steelblue","red"))+guides(fill="none")








