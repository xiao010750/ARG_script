#compostion
physeq_data <- physeq_arg  %>% 
  psmelt()%>%
  group_by(Type,Location,Sample) %>%
  summarise(total_count = sum(Abundance)) %>% #会创造一个新数据框，只包含分组变量和统计值。
  ungroup()#%>% 
#arrange(desc(total_count))

#t.test(total_count ~ Soil_type, data = physeq_data)

percent <- physeq_arg %>%
  psmelt()%>%
  group_by(Type) %>%
  summarise(total_count = sum(Abundance)) %>%
  ungroup()%>%
  group_by(Type) %>%
  summarise(total_count = sum(total_count)) %>%
  mutate(percentage = total_count / sum(total_count) * 100) %>% 
  data.frame()%>% 
  arrange(desc(percentage))

high_abundance_Type <- percent[percent$percentage>2.5,]$Type

physeq_data$Type[!(physeq_data$Type %in% high_abundance_Type)] <- "others"

physeq_data$Location <- factor(physeq_data$Location, levels = c("ZJ", "JS", "GD", "GX", "CQ", "SC", "YN"))

col_phylum <- c( "#808080","#FFA500","#CD1076","#F5BE8F","#A020F0","#87CEFA","#FFC0CB","#EE6A50","#4F94CD","#EE7600","#0395c6","#4876FF")

physeq_data <- physeq_data %>% 
  mutate(Type = fct_reorder(Type, total_count))

physeq_data %>% ggplot()+
  geom_col(aes(x=Sample,y=total_count,fill=Type),position = 'stack')+
  scale_fill_manual(values=col_phylum)+
  theme_minimal()+
  scale_y_continuous(expand = c(0,0))+
  labs(y="Abundance (copy per 16S rRNA gene)",x="")+
  theme(panel.grid = element_blank() , 
        panel.background = element_rect(fill = 'white') , 
        panel.border = element_rect(fill = NA , color = "black" ,  size = 0.8 ,  linetype = "solid"))+
  scale_y_continuous(limits=c(0,0.3))+
  facet_wrap(.~Location ,  
             nrow = 1,scales = "free_x")
ggsave("ARG_composition.pdf",width =10,height = 6)


#relative compostion
p1 <- physeq_data %>% subset(Location%in%c("ZJ")) %>% 
  ggplot(aes(x="",y=total_count,fill=Type))+
  geom_col() +
  scale_fill_manual(values = col_arg,name="Type")+
  coord_polar(theta = "y")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none");p1
p1+p2+p3+p4+p5+p6+p7+plot_layout(nrow = 1)

ggsave("ARG_compostion_rel.pdf",width =10,height = 6)



#alpha diversity
library(agricolae)
source("C:/Users/xiao0/Desktop/Rscript/alpha_diversity.R")
ps_arg<-physeq_arg 
ps_arg@sam_data$Location <- factor(ps_arg@sam_data[[1]], levels = c("ZJ", "JS", "GD", "GX", "CQ", "SC", "YN"))
physeq <- ps_arg %>% binary_otu
t_otu<-as.data.frame(t(as(otu_table(physeq),"matrix")))
group<-meta(physeq)
phy <- NULL
alpha_group<-cbind(jintao::alpha(t_otu,phy),group[rownames(t_otu),])
alpha_group_arg <- alpha_group
alpha_group <- left_join(alpha_group,read_csv("理化性质.csv"),by=c("sample_name"))

richness_out <- anova_sig(alpha_group , alpha_group$Richness , alpha_group$Location)
richness_out$type <- "Richness" 


alpha_out <- rbind( richness_out)%>%rename_with(~"marker" , 2)%>% rename_with(~"Location" , 3)

df_long <- alpha_group %>% select(Richness,Location,Soil_type) %>% 
  pivot_longer( cols = Richness , names_to = "type" , values_to = "alpha_index")

df_long_all <- left_join(df_long , alpha_out , by = c("type" , "Location"))


df_long_all$Soil_type <- factor(df_long_all$Soil_type, levels = c( "Bulk","Rhizosphere"))
df_long_all$Location <- factor(df_long_all$Location, levels = c("ZJ", "JS", "GD", "GX", "CQ", "SC", "YN"))
ggplot(df_long_all , aes(x = interaction(Soil_type,Location),y=alpha_index,shape =Soil_type))+
  scale_fill_manual(values = col_Location)+
  geom_boxplot(aes(fill = Location))+
  geom_jitter(aes(x = interaction(Soil_type,Location)  , alpha_index),size = 1.4)+
  scale_shape_manual(values = c(17,16))+
  geom_text(aes(x = interaction(Soil_type,Location)  , y = max+sd , label = marker) , size = 4 , position =  position_dodge(0.6))+
  theme(panel.grid = element_blank() , 
        panel.background = element_rect(fill = 'white') , 
        panel.border = element_rect(fill = NA , color = "black" ,  size = 0.8 ,  linetype = "solid"),
        axis.title = element_blank(), 
        axis.text.x = element_blank(),
        legend.position="none")

#beta
library(ggalt)
ps_arg<-physeq_arg
ps_arg@sam_data$Location <- factor(ps_arg@sam_data[[1]], levels = c("ZJ", "JS", "GD", "GX", "CQ", "SC", "YN"))
group<-meta(ps_arg)
otu<-otu_table(ps_arg) %>% data.frame()
distance <- vegdist(t(otu), method = 'bray')
nmds <- metaMDS(distance, k = 2)
# stress <- nmds$stress
adonis <- adonis2(t(otu)~Location+Soil_type,group);adonis

pp<-plot_ordination(ps_arg, ordinate(ps_arg, method = "NMDS", distance = "bray") , type="samples",color="Location",shape ="Soil_type")+geom_point(size=3)+
  scale_color_manual(values = col_Location)+
  scale_shape_manual(values = c(17,16))
pp+theme_bw()+#stat_ellipse()
  scale_x_continuous(expand = expand_scale(0.05))+scale_y_continuous(expand = expand_scale(0.05))+#expand_scale(0.05)表示将轴的数据范围扩展5%，这样可以留出一定空白。
  geom_encircle(aes(group =               Location,fill=Location),expand=0,spread=0.5,s_shape=1,size=3,linetype = 1,alpha=0.2)+
  scale_fill_manual(values = col_Location)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="bottom")
ggsave("ARG_beta.pdf",width = 6,height = 4)



# arg risk
col_arg <- c("#ef1828","#f88421","#ffbc14","#808080")
ARG_rank<-read.delim("E:/metagenomic/metagenomic_analysis/unique_gene/Sample_ranking_results.txt")
factor <- read.xlsx("group_factor.xlsx", sheetIndex = 1)
group_meta<-read.csv("group_metagenome.csv",row.names = 1) 
rank <- ARG_rank[,c(1:5)]
rank1 <- rank %>%
  mutate(
    Rank_I = Rank_I_per / (Rank_I_per + Rank_II_per + Rank_III_per + Rank_IV_per),
    Rank_II = Rank_II_per / (Rank_I_per + Rank_II_per + Rank_III_per + Rank_IV_per))

rank2 <- rank1 %>% select(Sample,Rank_I,Rank_II) %>% 
  pivot_longer( cols = -Sample , names_to = "type" , values_to = "Percent")
rank2$Sample <- gsub("\\.R_1.*", "", rank2$Sample)
colnames(group_meta)[3] <- c("Sample")
rank2 <- left_join(rank2,group_meta,by=c("Sample"))
rank2$Location <- factor(rank2$Location,levels = c("ZJ","JS","GD","GX","CQ","SC","YN"))
rank2$Sample <- factor(rank2$Sample,levels=factor$sample_name)
rank2 %>% 
  ggplot(aes(x=Sample,y=Percent,fill=type))+
  geom_col()+
  scale_fill_manual(values = col_arg)+
  scale_y_continuous(expand = c(0,0))+
  labs(y="Value")+
  theme_minimal()+
  theme(panel.grid = element_blank() , 
        panel.background = element_rect(fill = 'white') , 
        panel.border = element_rect(fill = NA , color = "black" ,  size = 0.8 ,  linetype = "solid"),
        axis.text.x = element_text(angle=45,hjust = 0.8,size = 7))+
  facet_wrap(.~Location ,  
             nrow = 1,scales = "free_x")
ggsave("ARG_risk.pdf",width =10,height = 6)

#venn plot
library(ggVennDiagram)
# rank="HMM.category" 
# physeq <- physeq_arg %>% speedyseq::tax_glom(taxrank = rank)
physeq <-physeq_arg %>%
  #subset_samples(Soil_type%in% c("Bulk")) %>% 
  binary_otu() %>% merge_samples2("Location",sum,unique)%>% 
  prune_taxa(taxa_sums(.)>0,.)

otu<- otu_table(physeq) %>% data.frame();otu[otu!=0]<-1 
res<-otu  %>% rownames_to_column("contigs") %>% melt() %>% filter(value!=0) 
res<-split(res$contigs,res$variable)

ordered_res <- res[c("ZJ", "JS", "GD", "GX", "CQ", "SC", "YN")]

a <- psmelt(physeq)
b <- a[,c(3,4,9)] %>% filter(Abundance!=0)
b<-split(b$HMM.category,b$Location)
ordered_b <- b[c("ZJ", "JS", "GD", "GX", "CQ", "SC", "YN")]
ggVennDiagram(ordered_b,label = "count",label_alpha = 0) +
  scale_fill_gradient(low = "#f8f8ff",high = "#f8f8ff")
ggsave("arg_venn.pdf",width = 6,height = 6)

#NST
library(NST)
col_com<- c("#ee9191","#4F94CD")
physeq <-physeq_arg
otu <- physeq@otu_table %>% data.frame() %>% t()
group <- meta(physeq)
group <- group[, c(1, 3)]
tnst=tNST(comm=otu, group=group, dist.method="bray",
          abundance.weighted=TRUE, rand=999,
          nworker=16, null.model="PF", between.group=TRUE,
          SES=TRUE, RC=TRUE)

result_nst<- tnst[[3]] %>% 
  group_by(group) %>% 
  summarise(value=mean(NST.ij.bray)) %>% 
  data.frame()
result_nst <-result_nst %>% 
  mutate(n_value=1-value) %>% 
  pivot_longer(., cols = c(value,n_value), names_to = "total", values_to = "Value")

result_nst$group <- factor(result_nst$group, levels = c("ZJ", "JS", "GD", "GX", "CQ", "SC", "YN"))
result_nst %>% ggplot()+
  geom_col(aes(x=group,y=Value,fill=total))+
  scale_fill_manual(values = col_com)+
  theme_minimal()+
  theme(panel.grid = element_blank() , 
        panel.background = element_rect(fill = 'white') , 
        panel.border = element_rect(fill = NA , color = "black" ,  size = 0.8 ,  linetype = "solid"))+
  labs(y = "Percent", x = "") +
  scale_y_continuous(expand = c(0,0))
ggsave("arg_community.pdf",width = 6,height = 4)

#
env <- read.csv("理化性质.csv")

mge <-ps_mge %>% 
  psmelt() %>% 
  group_by(Location,sample_name) %>% 
  summarise(total_count = sum(Abundance)) %>% #会创造一个新数据框，只包含分组变量和统计值。
  ungroup()
#t.test(total_count ~ Soil_type, data = mge)

arg <- physeq_arg %>% 
  psmelt()%>%
  group_by(Location,sample_name) %>%
  summarise(total_count = sum(Abundance)) %>% 
  ungroup()

mg <- left_join(mge,env,by="sample_name")
lm_model1 <- lm(total_count ~ Pre23, data = mg)
lm_model1 %>% summary()


ggplot(mg, aes(x = Pre23, y = total_count,shape=compartment))+
  geom_point(alpha=0.5,color="#4F94CD")+
  scale_shape_manual(values = c(17,16))+
  geom_smooth(aes(group = 1),method="lm")+
  theme_bw()+
  labs(y="Abundance",x="precipitation")#+stat_cor()

ggsave("arg_Richness_p.pdf",width = 6,height = 4)

#differential analysis
library(DAtest)
physeq <- ps_mge 
physeq <- physeq %>% speedyseq::tax_glom(taxrank = "Subtype")
ig <- DA.per(physeq,"Soil_type",relative = F)
dim(ig)
ig_filter <- ig[ig$pval.adj < 0.05, ]
ordering_counts <- table(ig_filter$ordering)
ordering_counts 


ig_filter[ig_filter$ordering=="Bulk>Rhizosphere",]$Name %>%
  unique() %>% 
  write.csv(file="bulk_rich_mge.csv")

#Heatmap
col_soil_type <- c("Bulk" = "#b18c3e", "Rhizosphere"="#a6d751")
library(ComplexHeatmap)
library(circlize)
physeq <- ps_mge
group<-meta(physeq)
rank="Subtype"
rich_mge <- read.csv("New_rhi_rich_mge.csv")

physeq<-ps_mge %>%subset_taxa(Name%in%rich_mge$Name) %>% 
  select_tax_table(rank) %>% 
  speedyseq::tax_glom(taxrank = rank) 

otu<-otu_table(physeq) %>% as.matrix();rownames(otu)<-data.frame(tax_table(physeq))[,rank];  otu<-log10(otu); otu[is.infinite(otu)]<-NA
tax<-tax_table(physeq) %>% data.frame %>% mutate(Subtype=ifelse(str_detect(Subtype," "), str_sub(Subtype, str_locate_my_first(Genus, " ",F)) ,Subtype)) 

soil_type <- sample_data(physeq)$Soil_type

sorted_col_names <- group_meta$sample_name

otu <- otu[, sorted_col_names] 


Heatmap(otu,name = "Log-transformed Abundance",
        top_annotation = HeatmapAnnotation(Soil_type = group[colnames(otu),]$Soil_type,col = list(Soil_type=col_soil_type),show_legend = FALSE), #是否显示注释图例
        row_labels  = tax$Subtype,
        heatmap_legend_param=list(title_position = "lefttop-rot",legend_height = unit(6, "cm")),
        col=c("#daeeff","#b3daff","#8dc7ff","#4ca6ff","#0081ff","#004c9a"),
        na_col =  "#daeeff",
        cluster_columns = F,cluster_rows = F,clustering_distance_rows="pearson")






