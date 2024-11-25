################################################################################
#Quality Control 
################################################################################
#Quality control to be done on 4 tse: microbiota, active, genes and functions 
#Bacteriota is a derivative of microbiota 
#tse_enzyme is derived from tse_pathway 
################################################################################
# Load the global variable 
mae <- readRDS("mae.rds")
tse = mae[[1]]
tse_bacteria = mae[[2]]
tse_active = mae[[3]]
tse_pathway = mae[[4]]
tse_enzyme = mae[[5]]
tse_gene = mae[[6]]


################################################################################
tse <- addPerCellQC(tse)

library_size_season = plotColData(tse,"sum","Season", colour_by = "Season", point_size = 5, 
                                  shape_by ="Season") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#png(filename="figures/library_size_season_kraken2.png" ,units = 'in',width=9, height=6, res=1000)
print(library_size_season)
#dev.off()
################################################################################
################################################################################
##Library size distribution

p1 <- ggplot(as.data.frame(colData(tse))) +
  geom_histogram(aes(x = sum), color = "black", fill = "gray", bins = 30) +
  labs(x = "Library size", y = "Frequency (n)") + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis


df_quality <- as.data.frame(colData(tse)) %>%
  arrange(sum) %>%
  mutate(index = 1:n())
p2 <- ggplot(df_quality, aes(y = index, x = sum/1e6)) +
  geom_point() +  
  labs(x = "Library size (million reads)", y = "Sample index") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis

#The distribution of calculated library sizes can be visualized as a histogram (left)
#or by sorting the samples by library size (right).

LibrarySize =p1 + p2

#png(filename="figures/library_size_kraken2.png" ,units = 'in',width=9, height=6, res=1000)
#print(LibrarySize)
#dev.off() 
################################
#Quality Control SAMSA2 
tse_active <- addPerCellQC(tse_active)

library_size_season_active = plotColData(tse_active,"sum","Season", colour_by = "Season", point_size = 5, 
                                         shape_by ="Season") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#png(filename="figures/library_size_season_active.png" ,units = 'in',width=9, height=6, res=1000)
#print(library_size_season_active)
#dev.off()
################################################################################
##Library size distribution

p1 <- ggplot(as.data.frame(colData(tse_active))) +
  geom_histogram(aes(x = sum), color = "black", fill = "gray", bins = 30) +
  labs(x = "Library size", y = "Frequency (n)") + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis


df_quality <- as.data.frame(colData(tse_active)) %>%
  arrange(sum) %>%
  mutate(index = 1:n())
p2 <- ggplot(df_quality, aes(y = index, x = sum/1e6)) +
  geom_point() +  
  labs(x = "Library size (million reads)", y = "Sample index") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis

#The distribution of calculated library sizes can be visualized as a histogram (left)
#or by sorting the samples by library size (right).

LibrarySize_Active =p1 + p2

#png(filename="figures/library_size_samsa2.png" ,units = 'in',width=9, height=6, res=1000)
#print(LibrarySize_Active)
#dev.off() 

##################################################
#Quality Control SAMSA2 functions
tse_pathway <- addPerCellQC(tse_pathway)

library_size_season_pathway = plotColData(tse_pathway,"sum","Season", colour_by = "Season", point_size = 5, 
                                          shape_by ="Season") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#png(filename="figures/library_size_season_pathway.png" ,units = 'in',width=9, height=6, res=1000)
#print(library_size_season_pathway)
#dev.off()
################################################################################
##Library size distribution

p1 <- ggplot(as.data.frame(colData(tse_pathway))) +
  geom_histogram(aes(x = sum), color = "black", fill = "gray", bins = 30) +
  labs(x = "Library size", y = "Frequency (n)") + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis


df_quality <- as.data.frame(colData(tse_pathway)) %>%
  arrange(sum) %>%
  mutate(index = 1:n())
p2 <- ggplot(df_quality, aes(y = index, x = sum/1e6)) +
  geom_point() +  
  labs(x = "Library size (million reads)", y = "Sample index") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), # Removes the grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) # Adds y-axis

#The distribution of calculated library sizes can be visualized as a histogram (left)
#or by sorting the samples by library size (right).

LibrarySize_pathway =p1 + p2

#png(filename="figures/library_size_pathway.png" ,units = 'in',width=9, height=6, res=1000)
#print(LibrarySize_pathway)
#dev.off() 

##################################################

