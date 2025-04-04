# Figure 1
pacman::p_load(
  ggplot2,
  dplyr,
  sf,
  arcgislayers
)

# Region color palette ----
getPalette <- c("#f4766c", "#ffd712",
                "#ef9b95", "#0a93d3",
                "#e50767", "#71e244",
                "#008121", "#ff80d9",
                "#858585", "#8999ff",
                "#0000b0", "#a10052",
                "#922d92", "#c87dff",
                "#ff8a00", "#00da94")
names(getPalette) <- c("CENTRAL-CORRIDOR-AF", "CENTRAL-AF",
                       "CENTRAL-PK", "EAST-PK",
                       "CENTRAL-CORRIDOR-PK", "GB",
                       "KARACHI", "KP",
                       "NORTH-AF", "NORTH-CORRIDOR-AF",
                       "NORTH-CORRIDOR-PK", "SINDH",
                       "SOUTH-CORRIDOR-AF", "SOUTH-CORRIDOR-PK",
                       "SOUTH-PUNJAB", "WEST-AF")


# A ----
# Annual stripes
low_season <- matrix(c("2006-01-01", "2007-01-01",
                       "2008-01-01", "2009-01-01",
                       "2010-01-01", "2011-01-01",
                       "2012-01-01", "2013-01-01",
                       "2014-01-01", "2015-01-01",
                       "2016-01-01", "2017-01-01",
                       "2018-01-01", "2019-01-01",
                       "2020-01-01", "2021-01-01",
                       "2022-01-01", "2023-01-01"),ncol=2,byrow=TRUE)

colnames(low_season) <- c("Start","End")
low_season %<>% as_tibble() %>%
  mutate(across(everything(), as.Date))

load("./Plot_data/Figure1A.RData") # pre-summarised data - case data is extracted from POLIS and linked to region

ggplot(WT1_AFP_region)+
  geom_col(aes(x = Y3M, y = cases, fill = epiblock), width = 92, show.legend = F) +
  geom_rect(data = low_season, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.3) +
  scale_fill_manual(values = getPalette) +
  coord_cartesian(xlim=  c(as.Date(date_decimal(2005.9)), as.Date(date_decimal(2023.7+.2))), ylim = c(0,130), clip="off")+
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Date", y = "Reported Cases")+
  theme_classic()
  
# B ----
# Map is plotted with the WHO polio administrative boundaries shapefile
# URL for WHO polio shapefiles
furl <- "https://services.arcgis.com/5T5nSi527N4F7luB/arcgis/rest/services/POLIO_ADMINISTRATIVE_BOUNDARIES/FeatureServer/4" 
# layers: 
# 0 = disputed borders 
# 1 = Disputed areas 
# 2 = ADM4 
# 3 = ADM3 
# 4 = ADM2 
# 5 = ADM1 
# 6 = ADM0 
# Add these to the end of the URL 

# Open connection
country_fl <- arc_open(furl)
# Subset with SQL
adm2_Pk_Af <- arc_select(country_fl,
                         where = "ISO_2_CODE = 'AF' or ISO_2_CODE = 'PK'")

# Only current regions
adm2_Pk_Af <- adm2_Pk_Af[as.Date(adm2_Pk_Af$ENDDATE)>=Sys.Date(),]

## Blocks assigned
adm2_Pk_Af$epiblock <- adm2_Pk_Af$ADM2_NAME

# Central corridor - Afghanistan
adm2_Pk_Af[adm2_Pk_Af$ADM1_NAME == "KHOST"
           | adm2_Pk_Af$ADM1_NAME == "PAKTIKA"
           | adm2_Pk_Af$ADM1_NAME == "PAKTYA" 
           | adm2_Pk_Af$ADM1_NAME == "GHAZNI" 
           ,]$epiblock <- "CENTRAL-CORRIDOR-AF"

# Central corridor - Pakistan
adm2_Pk_Af[adm2_Pk_Af$ADM2_NAME == "WAZIR-N" # KPTD
           | adm2_Pk_Af$ADM2_NAME == "WAZIR-S"
           | adm2_Pk_Af$ADM2_NAME == "FR BANNU" 
           | adm2_Pk_Af$ADM2_NAME == "FR DIKHAN"
           | adm2_Pk_Af$ADM2_NAME == "FR LAKKI"
           | adm2_Pk_Af$ADM2_NAME == "FR TANK"
           
           | adm2_Pk_Af$ADM2_NAME == "BANNU"      # KP
           | adm2_Pk_Af$ADM2_NAME == "DIKHAN"
           
           | adm2_Pk_Af$ADM2_NAME == "LAKKIMRWT" 
           | adm2_Pk_Af$ADM2_NAME == "TANK"
           ,]$epiblock <- "CENTRAL-CORRIDOR-PK"


# centre_pk - Balochistan
adm2_Pk_Af[adm2_Pk_Af$ADM2_NAME == "BARKHAN"       
           | adm2_Pk_Af$ADM2_NAME == "BOLAN"
           | adm2_Pk_Af$ADM2_NAME == "DBUGTI" 
           | adm2_Pk_Af$ADM2_NAME == "DUKKI"        # Loralai
           | adm2_Pk_Af$ADM2_NAME == "HARNAI"
           | adm2_Pk_Af$ADM2_NAME == "JAFARABAD"
           | adm2_Pk_Af$ADM2_NAME == "JHALMAGSI"
           | adm2_Pk_Af$ADM2_NAME == "KALAT"
           | adm2_Pk_Af$ADM2_NAME == "SSIKANDARABAD" # KALAT
           | adm2_Pk_Af$ADM2_NAME == "KHUZDAR"
           | adm2_Pk_Af$ADM2_NAME == "KOHLU"
           | adm2_Pk_Af$ADM2_NAME == "KSAIFULAH"  
           | adm2_Pk_Af$ADM2_NAME == "LASBELA"
           | adm2_Pk_Af$ADM2_NAME == "LORALAI" 
           | adm2_Pk_Af$ADM2_NAME == "MASTUNG"
           | adm2_Pk_Af$ADM2_NAME == "NSIRABAD" 
           | adm2_Pk_Af$ADM2_NAME == "SHARANI"
           | adm2_Pk_Af$ADM2_NAME == "SIBI"
           | adm2_Pk_Af$ADM2_NAME == "SOHBATPUR"    # JAFARABAD
           | adm2_Pk_Af$ADM2_NAME == "ZHOB"  
           | adm2_Pk_Af$ADM2_NAME == "ZIARAT"
           ,]$epiblock <- "CENTRAL-PK"

# province name in both countries
adm2_Pk_Af[adm2_Pk_Af$ADM1_NAME == "BALOCHISTAN" & adm2_Pk_Af$ADM2_NAME == "MUSAKHEL",]$epiblock <-  "CENTRAL-PK"

# Sindh - full Sindh
adm2_Pk_Af[adm2_Pk_Af$ADM1_NAME == "SINDH"       # Remove Karachi later
           ,]$epiblock <- "SINDH"
# east_af
adm2_Pk_Af[adm2_Pk_Af$ADM1_NAME == "BAMYAN"
           | adm2_Pk_Af$ADM1_NAME == "KABUL" 
           | adm2_Pk_Af$ADM1_NAME == "KAPISA" 
           | adm2_Pk_Af$ADM1_NAME == "LOGAR" 
           | adm2_Pk_Af$ADM1_NAME == "PANJSHER"
           | adm2_Pk_Af$ADM1_NAME == "PARWAN" 
           | adm2_Pk_Af$ADM1_NAME == "WARDAK" 
           | adm2_Pk_Af$ADM1_NAME == "DAYKUNDI" 
           ,]$epiblock <- "CENTRAL-AF"

# east_pk - full provinces
adm2_Pk_Af[adm2_Pk_Af$ADM1_NAME == "AJK" 
           | adm2_Pk_Af$ADM1_NAME == "GBALTISTAN" 
           | adm2_Pk_Af$ADM1_NAME == "ISLAMABAD" 
           | adm2_Pk_Af$ADM1_NAME == "PUNJAB"
           ,]$epiblock <- "EAST-PK"

adm2_Pk_Af[adm2_Pk_Af$ADM1_NAME == "GILGIT BALTISTAN",]$epiblock <- "GB"


#KP remainder
adm2_Pk_Af[adm2_Pk_Af$ADM2_NAME == "ABOTABAD"
           | adm2_Pk_Af$ADM2_NAME == "BATAGRAM"
           | adm2_Pk_Af$ADM2_NAME == "BUNER"
           | adm2_Pk_Af$ADM2_NAME == "HARIPUR"
           | adm2_Pk_Af$ADM2_NAME == "KOHISTAN LOWER"
           | adm2_Pk_Af$ADM2_NAME == "KOHISTAN UPPER"
           | adm2_Pk_Af$ADM2_NAME == "MANSEHRA"
           | adm2_Pk_Af$ADM2_NAME == "SHANGLA"
           | adm2_Pk_Af$ADM2_NAME == "TORGHAR"
           | adm2_Pk_Af$ADM2_NAME == "KOLAI PALAS"
           ,]$epiblock <- "KP"


# South Punjab
adm2_Pk_Af[adm2_Pk_Af$ADM2_NAME == "MULTAN"
           | adm2_Pk_Af$ADM2_NAME == "KHANEWAL"
           | adm2_Pk_Af$ADM2_NAME == "LODHRAN"
           | adm2_Pk_Af$ADM2_NAME == "VEHARI"
           | adm2_Pk_Af$ADM2_NAME == "BAHAWALPUR"
           | adm2_Pk_Af$ADM2_NAME == "RYKHAN"
           | adm2_Pk_Af$ADM2_NAME == "LAYYAH"
           | adm2_Pk_Af$ADM2_NAME == "RAJANPUR"
           | adm2_Pk_Af$ADM2_NAME == "PAKPATTEN"
           | adm2_Pk_Af$ADM2_NAME == "KHANEWAL"
           | adm2_Pk_Af$ADM2_NAME == "MUZFARGARH"
           | adm2_Pk_Af$ADM2_NAME == "BAHWLNAGAR"
           | adm2_Pk_Af$ADM2_NAME == "DGKHAN"
           
           ,]$epiblock <- "SOUTH-PUNJAB"


# north_af
adm2_Pk_Af[adm2_Pk_Af$ADM1_NAME == "BADAKHSHAN"
           | adm2_Pk_Af$ADM1_NAME == "BAGHLAN"
           | adm2_Pk_Af$ADM1_NAME == "BALKH" 
           | adm2_Pk_Af$ADM1_NAME == "FARYAB" 
           | adm2_Pk_Af$ADM1_NAME == "JAWZJAN" 
           | adm2_Pk_Af$ADM1_NAME == "KUNDUZ" 
           | adm2_Pk_Af$ADM1_NAME == "SAMANGAN"
           | adm2_Pk_Af$ADM1_NAME == "SARI PUL"    # both are present
           | adm2_Pk_Af$ADM1_NAME == "SAR-E-PUL"   
           | adm2_Pk_Af$ADM1_NAME == "TAKHAR" 
           ,]$epiblock <- "NORTH-AF"

# north_corridor - Afghanistan
adm2_Pk_Af[adm2_Pk_Af$ADM1_NAME == "KUNAR"
           | adm2_Pk_Af$ADM1_NAME == "NANGARHAR"
           | adm2_Pk_Af$ADM1_NAME == "NURISTAN" 
           | adm2_Pk_Af$ADM1_NAME == "LAGHMAN" 
           ,]$epiblock <- "NORTH-CORRIDOR-AF"

# north_corridor - Pakistan
adm2_Pk_Af[adm2_Pk_Af$ADM2_NAME == "BAJOUR"    # KPTD
           | adm2_Pk_Af$ADM2_NAME == "KHYBER"
           | adm2_Pk_Af$ADM2_NAME == "MOHMAND"
           | adm2_Pk_Af$ADM2_NAME == "FR PESHAWAR"
           | adm2_Pk_Af$ADM2_NAME == "ORAKZAI"
           | adm2_Pk_Af$ADM2_NAME == "KURRAM"
           | adm2_Pk_Af$ADM2_NAME == "FR KOHAT" 
           
           | adm2_Pk_Af$ADM2_NAME == "CHARSADA"  # KP
           | adm2_Pk_Af$ADM2_NAME == "CHITRAL" # Sometimes split to lower and upper
           | adm2_Pk_Af$ADM2_NAME == "CHITRAL LOWER"
           | adm2_Pk_Af$ADM2_NAME == "CHITRAL UPPER"
           | adm2_Pk_Af$ADM2_NAME == "DIRLOWER"
           | adm2_Pk_Af$ADM2_NAME == "DIRUPPER"
           | adm2_Pk_Af$ADM2_NAME == "MALAKAND"
           | adm2_Pk_Af$ADM2_NAME == "NOWSHERA"
           | adm2_Pk_Af$ADM2_NAME == "PESHAWAR"
           | adm2_Pk_Af$ADM2_NAME == "SWABI"
           | adm2_Pk_Af$ADM2_NAME == "SWAT"
           | adm2_Pk_Af$ADM2_NAME == "MARDAN"
           | adm2_Pk_Af$ADM2_NAME == "KOHAT" 
           | adm2_Pk_Af$ADM2_NAME == "HANGU" 
           | adm2_Pk_Af$ADM2_NAME == "KARAK"
           ,]$epiblock <- "NORTH-CORRIDOR-PK"

# south_corridor - Afghanistan
adm2_Pk_Af[adm2_Pk_Af$ADM1_NAME == "HILMAND"
           | adm2_Pk_Af$ADM1_NAME == "KANDAHAR"
           | adm2_Pk_Af$ADM1_NAME == "NIMROZ"
           | adm2_Pk_Af$ADM1_NAME == "URUZGAN"
           | adm2_Pk_Af$ADM1_NAME == "ZABUL" 
           ,]$epiblock <- "SOUTH-CORRIDOR-AF"

# south_corridor - Pakistan
adm2_Pk_Af[adm2_Pk_Af$ADM2_NAME == "KABDULAH"    
           | adm2_Pk_Af$ADM2_NAME == "PISHIN"
           | adm2_Pk_Af$ADM2_NAME == "QUETTA"
           ,]$epiblock <- "SOUTH-CORRIDOR-PK"

# west_af
adm2_Pk_Af[ adm2_Pk_Af$ADM1_NAME == "BADGHIS"
            | adm2_Pk_Af$ADM1_NAME == "FARAH"
            | adm2_Pk_Af$ADM1_NAME == "GHOR"
            | adm2_Pk_Af$ADM1_NAME == "HIRAT"
            ,]$epiblock <- "WEST-AF"

# west_pk
adm2_Pk_Af[adm2_Pk_Af$ADM2_NAME == "AWARAN"    
           | adm2_Pk_Af$ADM2_NAME == "CHAGHAI"
           | adm2_Pk_Af$ADM2_NAME == "GWADUR"
           | adm2_Pk_Af$ADM2_NAME == "KECH"
           | adm2_Pk_Af$ADM2_NAME == "KHARAN"
           | adm2_Pk_Af$ADM2_NAME == "NOSHKI"
           | adm2_Pk_Af$ADM2_NAME == "PANJGOUR"
           | adm2_Pk_Af$ADM2_NAME == "WASHUK"
           ,]$epiblock <- "CENTRAL-PK"

# Karachi
adm2_Pk_Af[adm2_Pk_Af$ISO_2_CODE == "PK" & grepl("KHI", adm2_Pk_Af$ADM2_NAME),]$epiblock <- "KARACHI"


## Plot
## plot with country borders
furl_0 <- "https://services.arcgis.com/5T5nSi527N4F7luB/arcgis/rest/services/POLIO_ADMINISTRATIVE_BOUNDARIES/FeatureServer/6" 
# Open connection
country_fl_0 <- arc_open(furl_0)
# Subset with SQL
adm0_Pk_Af <- arc_select(country_fl_0,
                         where = "ISO_2_CODE = 'AF' or ISO_2_CODE = 'PK'")


adm0_Pk_Af <- adm0_Pk_Af[as.Date(adm0_Pk_Af$ENDDATE)>=Sys.Date(),]

# Plot
ggplot() +
  geom_sf(data = adm2_Pk_Af, aes(fill = epiblock), color = "lightgrey", linewidth = 0.2)+
  geom_sf(data = adm0_Pk_Af, fill = NA, linewidth = 0.7, color = "black")+
  scale_fill_manual(values = getPalette)+
  theme_void()

# published plot has some additional aesthetic elements not recreated here for ease of use

