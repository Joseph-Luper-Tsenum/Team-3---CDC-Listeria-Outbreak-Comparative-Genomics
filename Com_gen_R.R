args <- commandArgs(trailingOnly = TRUE)
rnorm(MD=args[1], kSNPtree=args[2], CHEWtree=args[3], STRINGtree=args[4], ANItree=args[5])

#### Housekeeping ####
library(ape)
library(ggtree)
library(ggplot2)
library(usmap)
library(phylotools)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtree")
theme_alpha <- theme_bw() + theme(axis.title.x = element_text(size = 16),
                                  axis.title.y = element_text(size = 16),
                                  axis.text = element_text(size = 12),
                                  strip.text = element_text(size=16),
                                  panel.grid = element_blank(),
                                  legend.position = "right", legend.title = element_blank(),
                                  legend.text = element_text(size = 12))
#### Import Metadata ####
team3FoodConsumptionQ3689 <- read.delim(MD, header=FALSE)
team3FoodConsumptionQ3689$V8 <- NULL
nms <- c("ID", "Location", "SampleDate", "Food1", "Food2", "Food3", "Food4")
colnames(team3FoodConsumptionQ3689) <- nms

list <- unique(team3FoodConsumptionQ3689$Food1)
list <- append(list,unique(team3FoodConsumptionQ3689$Food2))
list <- append(list,unique(team3FoodConsumptionQ3689$Food3))
list <- append(list,unique(team3FoodConsumptionQ3689$Food4))

foods <- unique(list)
foods <- foods[-1]
length(foods)
dim(team3FoodConsumptionQ3689)

df <- data.frame(matrix(NA,nrow = 50,ncol = 24))
b <- c(colnames(team3FoodConsumptionQ3689)[1:3],foods)
colnames(df) <- b

df$ID <- team3FoodConsumptionQ3689$ID
df$Location <- team3FoodConsumptionQ3689$Location
df$SampleDate <- team3FoodConsumptionQ3689$SampleDate

food.grep <- paste0(foods, '$')

for(i in 1:length(food.grep)){
  for(j in 1:length(df$ID)){
    temp = sapply(team3FoodConsumptionQ3689[j,c(4:7)], function(x) grepl(food.grep[i], x))
    c <- length(temp[temp== TRUE])
    df[j,foods[i]] <- c
  }
}
MD <- df
MD$date_num <- as.POSIXct(MD$SampleDate, format="%Y-%m-%d")
MD$date_num <- as.Date(MD$date_num)
row.names(MD) <- MD$ID
MD$X <- NULL
MD$Queso.Fresco <- MD$Queso.Fresco + MD$Queso.Fresco.
MD$Queso.Fresco. <- NULL
for(i in 4:23){
  MD[,i] <- as.factor(MD[,i])
  levels(MD[,i]) <- c("Absent", "Present")
}

for(i in 1:length(MD$Location)){
  if(MD$Location[i]=="Virginia"){
    MD$latitude[i] = 37.926868
    MD$longitude[i] = -78.024902
  }else if(MD$Location[i]=="Kansas"){
    MD$latitude[i] = 38.500000
    MD$longitude[i] = -98.000000
  }else if(MD$Location[i]=="Maryland"){
    MD$latitude[i] = 39.045753
    MD$longitude[i] = -76.641273
  }else if(MD$Location[i]=="Alabama"){
    MD$latitude[i] = 32.318230
    MD$longitude[i] = -86.902298
  }else if(MD$Location[i]=="New York"){
    MD$latitude[i] = 40.730610
    MD$longitude[i] = -73.935242
  }else if(MD$Location[i]=="Ohio"){
    MD$latitude[i] = 40.367474
    MD$longitude[i] = -82.996216
  }else if(MD$Location[i]=="South Carolina"){
    MD$latitude[i] = 33.836082
    MD$longitude[i] = -81.163727
  }
}

#### Import trees####
tre.resolved <- read.tree(kSNPtree)
tre.fastANI.NJ <- read.tree(ANItree)
tre.chew.NJ <- read.tree(CHEWtree)

names <- data.frame(old.lab = tre.chew.NJ$tip.label ,new.lab = substr(tre.chew.NJ$tip.label,start =1,stop=7))
tre.chew.NJ <- phylotools::sub.taxa.label(tre.chew.NJ, names)

tre.stringmlst.NJ <- read.tree(STRINGtree)

#### Food tree ####
P <- ggtree(tre.resolved, layout = "rectangular") %<+% MD 
# suspicious:
# Queso.Fresco
Quesoplot <- P + geom_tiplab(hjust = -0.1)+ 
  geom_tippoint(aes(color = Queso.Fresco), size=5) + 
  theme(legend.position = "right") +  
  geom_treescale() +
  scale_color_manual(values=c("#999999", "darkred")) +
  labs(color = "Queso Fresco")
Quesoplot
ggsave("Quesoplot.jpeg",Quesoplot,width = 20, height = 10)
# Cream.Cheese
CreamCheeseplot <- P + geom_tiplab(hjust = -0.1)+ 
  geom_tippoint(aes(color = Cream.Cheese), size=5) + 
  theme(legend.position = "right") +  
  geom_treescale() +
  scale_color_manual(values=c("#999999", "darkred")) +
  labs(color = "Cream Cheese")
CreamCheeseplot
ggsave("CreamCheeseplot.jpeg",CreamCheeseplot,width = 20, height = 10)
# Greek.Salad
GreekSaladplot <- P + geom_tiplab(hjust = -0.1)+ 
  geom_tippoint(aes(color = Greek.Salad), size=5) + 
  theme(legend.position = "right") +  
  geom_treescale() +
  scale_color_manual(values=c("#999999", "darkred")) +
  labs(color = "Greek Salad")
GreekSaladplot
ggsave("GreekSaladplot.jpeg",GreekSaladplot,width = 20, height = 10)

#### Outbreak analysis ####

outbreak = c('CGT1239', 'CGT1729', 'CGT1144','CGT1624','CGT1743','CGT1782', 
             'CGT1288', 'CGT1700', 'CGT1632', 'CGT1113', 
             'CGT1627', 'CGT1540', 'CGT1029', 'CGT1987', 'CGT1121', 'CGT1521', 
             'CGT1756','CGT1376','CGT1832', 'CGT1687','CGT1111') # 'CGT1939','CGT1579',

for(i in 1:length(MD$ID)){
  if(MD$ID[i] %in% outbreak){
    MD$Outbreak[i] <- 'Outbreak'
  }else{
    MD$Outbreak[i] <- 'Sporadic'
  }
}

# visualize outbreak over time
p <- ggplot(MD, aes(x=date_num, fill=Outbreak)) +
  geom_histogram(bins = 40)+theme_alpha+xlab("Time")+scale_x_date(date_labels = "%Y %b %d")+
  scale_fill_manual(values=c("darkred","#999999"))
p 
ggsave("outbreak_over_time.jpeg",p,height = 7, width = 7)


# Make fastANI tree
P <- ggtree(tre.fastANI, layout = "rectangular") %<+% MD 
Outbreak.tree <- P + geom_tiplab(hjust = -0.1)+ 
  geom_tippoint(aes(color = Outbreak), size=5) + 
  theme(legend.position = "right") +  
  geom_treescale() +
  scale_color_manual(values=c("darkred", "#999999")) +
  labs(color = "Outbreak")
Outbreak.tree
ggsave("tree_fastANI.jpeg",Outbreak.tree,width = 20, height = 10)

# make the "resolved" SNP tree
P <- ggtree(tre.resolved, layout = "rectangular") %<+% MD 
Outbreak.tree <- P + geom_tiplab(hjust = -0.1)+ 
  geom_tippoint(aes(color = Outbreak), size=5) + 
  theme(legend.position = "right") +  
  geom_treescale() +
  scale_color_manual(values=c("darkred", "#999999")) +
  labs(color = "Outbreak")
Outbreak.tree

ggsave("tree_RESOLVED.jpeg",Outbreak.tree,width = 20, height = 10)

# make the sting MLST tree
P <- ggtree(tre.stringmlst.NJ, layout = "rectangular") %<+% MD 
Outbreak.tree <- P + geom_tiplab(hjust = -0.1)+ 
  geom_tippoint(aes(color = Outbreak), size=5) + 
  theme(legend.position = "right") +  
  geom_treescale() +
  scale_color_manual(values=c("darkred", "#999999")) +
  labs(color = "Outbreak")
Outbreak.tree

ggsave("tree_stingMLST.jpeg",Outbreak.tree,width = 20, height = 10)

# make the Chebacca tree
P <- ggtree(tre.chew.NJ, layout = "rectangular") %<+% MD 
Outbreak.tree <- P + geom_tiplab(hjust = -0.1)+ 
  geom_tippoint(aes(color = Outbreak), size=5) + 
  theme(legend.position = "right") +  
  geom_treescale() +
  scale_color_manual(values=c("darkred", "#999999")) +
  labs(color = "Outbreak")
Outbreak.tree

ggsave("tree_Chebacca.jpeg",Outbreak.tree,width = 20, height = 10)

#### MAPPING ####
loc <- data.frame(lon = MD$longitude, lat = MD$latitude )

transformed_data <- usmap_transform(loc)

for(i in 1:length(MD$Location)){
  if(MD$latitude[i]==37.926868){
    MD$x[i] = 1903574.9
    MD$y[i] = -532775.6
  }else if(MD$latitude[i]==38.500000){
    MD$x[i] = 174303.7
    MD$y[i] = -720288.9
  }else if(MD$latitude[i]==39.045753){
    MD$x[i] = 1987011.1
    MD$y[i] = -378939.7
  }else if(MD$latitude[i]==32.318230){
    MD$x[i] = 1232466.4
    MD$y[i] = -1312794.8
  }else if(MD$latitude[i]==40.730610){
    MD$x[i] = 2152333.9
    MD$y[i] = -128964.0
  }else if(MD$latitude[i]==40.367474){
    MD$x[i] = 1429122.4
    MD$y[i] = -366976.3
  }else if(MD$latitude[i]==33.836082){
    MD$x[i] = 1730524.6
    MD$y[i] = -1046391.5
  }
}

MD.filt <- MD[MD$Outbreak=="Outbreak",]
map <- plot_usmap("states") + 
  geom_count(data = MD, 
             aes(x = x, y = y), 
             color = "#999999")+
  geom_count(data = MD.filt, 
             aes(x = x, y = y), 
             color = "darkred")

ggsave("Outbreak_map.jpeg", map)
