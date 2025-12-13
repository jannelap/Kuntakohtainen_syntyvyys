library("dplyr")
library("tidyr")
library("sf")
library("geofi")
library("ggplot2")
library("rstan")
library("dynamite")
library("bayesplot")
library("devtools")
library("posterior")
library("latex2exp")
library("stringr")
library("tibble")
#
#------------------------------------------------
# Kuntakohtaisen syntyvyyden mallintaminen
# dynaamisella monimuuttujaisella paneelimallilla
#
# Tekijä: Janne Lappalainen
#------------------------------------------------
#
#ggplot teema
theme_set(theme_bw() + theme(panel.background = element_rect(fill = 'white', colour = 'black'),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.1),
      panel.grid.minor = element_line(colour = "grey90")))

#-------------------------------------------
# Funktiot
#-------------------------------------------

#Laskee pienimmät tehokkaat otoskoot ja suurimmat R-hatut
ess_r <- function(malli){
  d<-posterior::as_draws_rvars(malli$stanfit)
  
  min_ess_b <- NA
  min_ess_t <- NA
  max_rhat <- NA
  
  for (i in 1:length(d)){
    bulk <- ess_bulk(d[[i]])
    tail <- ess_tail(d[[i]])
    r <- rhat(d[[i]])
    if(is.na(min_ess_b) | min(bulk) < min_ess_b){
      min_ess_b <- min(bulk)
    }
    if(is.na(min_ess_t) | min(tail) < min_ess_t){
      min_ess_t <- min(tail)
    }
    if(is.na(max_rhat) | max(r) > max_rhat){
      max_rhat <- max(r)
    }
  }
  return(c(min_ess_b, min_ess_t, max_rhat))
}


#Laskee selitysasteen R2
R2 <- function(naytteet){
  naytteet_r2 <- naytteet %>%
    filter(Vuosi>1990) %>%
      mutate(Jaannos = Hedelmallisyys - Hedelmallisyys_fitted) %>%
        group_by(.draw) %>%
          summarize(Hedelmallisyys_fitted_var = var(Hedelmallisyys_fitted), Jaannos_var = var(Jaannos))
  
  naytteet_r2 <- naytteet_r2 %>%
    mutate(R2 = Hedelmallisyys_fitted_var/(Hedelmallisyys_fitted_var + Jaannos_var))
  
  return(naytteet_r2$R2)
}

#Laskee aikariippuvien parametrien posteriorinäytteiden keskiarvot, 95%:n posteriorivälit,
#positiivisten ja negatiivisten arvojen vuosittaiset todennäköisyydet,
#todennäköisyyden, että arvo on positiivinen joka vuonna ja todennäköisyyden, että arvo on negatiivinen joka vuonna
posteriori_tunnusluvut <- function(naytteet){
  #Käännetään posteriorinäytteiden matriisi
  naytteet<-as.data.frame(t(naytteet))
  
  #Otosten sarakkeiden nimet
  otokset <- colnames(naytteet)
  
  #Luodaan sarake selittäjälle ja vuodelle
  naytteet <- naytteet %>%
    rownames_to_column(var = "Selittaja") %>%
      filter(!Selittaja %in% c(".chain", ".iteration", ".draw")) %>%
        mutate(Vuosi = str_extract(Selittaja, "(19|20)\\d{2}"),
        Selittaja = str_remove(Selittaja, "\\[(19|20)\\d{2}\\]"))
  
  #Lasketaan posteriorikeskiarvo ja 95% posteriorivälit
  naytteet <- naytteet %>%
    rowwise() %>%
      mutate(Ka = mean(c_across(all_of(otokset))), Ala95 = quantile(c_across(all_of(otokset)), probs = 0.025), Yla95 = quantile(c_across(all_of(otokset)), probs = 0.975))
  
  #Lasketaan positiivisten näytteiden osuus
  naytteet <- naytteet %>%
    mutate(across(all_of(otokset), ~ ifelse(.>0, 1, 0))) %>%
      rowwise() %>%
        mutate(Positiiviset_ka = mean(c_across(all_of(otokset))))
  
  #Lasketaan osuus niille näytteille, joilla on jokaisena vuotena positiivinen arvo
  #ja niille, joilla jokaisena vuotena on negatiivinen arvo
  naytteet_summa <- naytteet %>%
    group_by(Selittaja) %>%
      summarise(across(all_of(otokset), mean)) %>%
          mutate(across(all_of(otokset), ~ ifelse(.==1, 1, 0), .names = "Posi_{.col}")) %>%
            mutate(across(all_of(otokset), ~ ifelse(.==0, 1, 0), .names = "Nega_{.col}")) %>%
              rowwise() %>%
                mutate(Kaikki_positiiviset = mean(c_across(starts_with("Posi_"))),
                Kaikki_negatiiviset = mean(c_across(starts_with("Nega_")))) %>%
                  select(Selittaja, Kaikki_positiiviset, Kaikki_negatiiviset)

  #Yhdistetään tulokset
  naytteet_tunnusluvut <- left_join(naytteet, naytteet_summa, by = "Selittaja") %>%
    mutate(Negatiiviset_ka = 1 - Positiiviset_ka) %>%
      select(Selittaja, Vuosi, Ka, Ala95, Yla95, Positiiviset_ka, Negatiiviset_ka, Kaikki_positiiviset, Kaikki_negatiiviset)
  
  return(naytteet_tunnusluvut)
}

#-------------------------------------------
# Datan käsittely
#-------------------------------------------

#Lasketaan kokonaishedelmällisyysluvut

#Elävänä syntyneet äidin iän mukaan
#12dq
nimet_syntyneet<-c("Kunta", "Vuosi", "lapset15_19", "lapset20_24", "lapset25_29", "lapset30_34","lapset35_39", "lapset40_44", "lapset45_49")
syntyneet<-read.csv2(file="ainestot/2025_aluejako/syntyneet.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_syntyneet )

#Naisten keskiväkiluku iän mukaan
#11s1
nimet_vakiluku<-c("Kunta", "Vuosi", "naiset15_19", "naiset20_24", "naiset25_29", "naiset30_34","naiset35_39", "naiset40_44", "naiset45_49")
vakiluku<-read.csv2(file="ainestot/2025_aluejako/keskivakiluku.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_vakiluku )

kunnat<-left_join(syntyneet, vakiluku, by = c("Kunta", "Vuosi"))

#Erotetaan kunnan tunnus ja nimi
kunnat <- kunnat %>%
  mutate(Id = as.integer(substring(Kunta, first = 1, last = 3)),
         Kunta = substring(Kunta, first = 5))

kunnat <- kunnat %>%
  mutate(naiset15_19 = as.double(naiset15_19),
         naiset20_24 = as.double(naiset20_24),
         naiset25_29 = as.double(naiset25_29),
         naiset30_34 = as.double(naiset30_34),
         naiset35_39 = as.double(naiset35_39),
         naiset40_44 = as.double(naiset40_44),
         naiset45_49 = as.double(naiset45_49)) 


#Lasketaan kuntien ikäluokkakohtaiset hedelmällisyysluvut
kunnat_afr <- kunnat %>%
  mutate(h15_19 = lapset15_19 / naiset15_19,
         h20_24 = lapset20_24 / naiset20_24,
         h25_29 = lapset25_29 / naiset25_29,
         h30_34 = lapset30_34 / naiset30_34,
         h35_39 = lapset35_39 / naiset35_39,
         h40_44 = lapset40_44 / naiset40_44,
         h45_49 = lapset45_49 / naiset45_49) 

#Lasketaan kuntien kokonaishedelmällisyysluvut
kunnat <- kunnat_afr %>%
  mutate(Hedelmallisyys = 5*(h15_19 + h20_24 + h25_29 + h30_34 + h35_39 + h40_44 + h45_49)) 

#Tarkistetaan puuttuvat hedelmällisyysluvut
puuttuvat <- kunnat[which(is.na(kunnat$Hedelmallisyys)),]
dim(puuttuvat)

#Pudotetaan Sottunga, koska siellä on ikäluokkia, joille ei saada laskettua hedelmällisyyslukua
kunnat<-kunnat %>%
  filter(Kunta != "Sottunga")


#Selittäjät
#Naisten osuus 19-49 vuotiaista
nimet_selittajat<-c("Kunta","Vuosi", "Yhteensa", "Naiset")
naiset_19_49<-read.csv2(file="ainestot/2025_aluejako/naiset_19_49.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_selittajat)

#Lasketaan suhteellinen osuus
naiset_19_49 <- naiset_19_49 %>%
  #Korjataan Maarianhaminan kirjoitusasu
  mutate(Kunta = if_else(Kunta == "Maarianhamina - Mariehamn", "Maarianhamina", Kunta),
         Naiset_osuus_19_49 = Naiset/Yhteensa
  ) %>%
  select(Kunta, Vuosi, Naiset_osuus_19_49)


#Väkiluku, Keski-ikä, Uskontokuntiin kuuumattomien osuus ja Asuinalueellaan syntyneiden osuus
nimet_asuinalue<-c("Kunta","Vuosi", "Vakiluku", "Keski_ika", "Ei_uskontokuntaa", "Asuinalueella_syntyneet")
tunnuslukuja<-read.csv2(file="ainestot/2025_aluejako/11ra.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_asuinalue)

tunnuslukuja <- tunnuslukuja %>%
  #Korjataan Maarianhaminan kirjoitusasu
  mutate(Kunta = if_else(Kunta == "Maarianhamina - Mariehamn", "Maarianhamina", Kunta),
         log_Vakiluku = log(Vakiluku)
  )


#Asuntokuntien koot
nimet_selittajat<-c("Kunta", "Vuosi","Yhteensa","Yhden_henkilon_asuntokunnat","Kahden_henkilon_asuntokunnat","Kolmen_henkilon_asuntokunnat","Neljan_henkilon_asuntokunnat","Viiden_henkilon_asuntokunnat","Kuuden_henkilon_asuntokunnat","Vahintaan_seitseman_henkilon_asuntokunnat")
asuntokunnat<-read.csv2(file="ainestot/2025_aluejako/asuntokuntien_koot.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_selittajat)

asuntokunnat <- asuntokunnat %>%
  #Korjataan Maarianhaminan kirjoitusasu
  mutate(Kunta = if_else(Kunta == "Maarianhamina - Mariehamn", "Maarianhamina", Kunta)
  )

asuntokunnat <- asuntokunnat %>%
  mutate(Yhden_henkilon_asuntokunnat_osuus = Yhden_henkilon_asuntokunnat/Yhteensa,
         Kahden_henkilon_asuntokunnat_osuus = Kahden_henkilon_asuntokunnat/Yhteensa,
         Kolmen_henkilon_asuntokunnat_osuus = Kolmen_henkilon_asuntokunnat/Yhteensa,
         Neljan_henkilon_asuntokunnat_osuus = Neljan_henkilon_asuntokunnat/Yhteensa,
         Viiden_henkilon_asuntokunnat_osuus = Viiden_henkilon_asuntokunnat/Yhteensa,
         Kuuden_henkilon_asuntokunnat_osuus = Kuuden_henkilon_asuntokunnat/Yhteensa,
         Vahintaan_kuuden_henkilon_asuntokunnat_osuus = (Kuuden_henkilon_asuntokunnat+Vahintaan_seitseman_henkilon_asuntokunnat)/Yhteensa,
         Vahintaan_seitseman_henkilon_asuntokunnat_osuus = Vahintaan_seitseman_henkilon_asuntokunnat/Yhteensa
  )%>%
  select(Kunta, Vuosi, Yhden_henkilon_asuntokunnat_osuus, Kahden_henkilon_asuntokunnat_osuus, Kolmen_henkilon_asuntokunnat_osuus, Neljan_henkilon_asuntokunnat_osuus, Viiden_henkilon_asuntokunnat_osuus, Kuuden_henkilon_asuntokunnat_osuus, Vahintaan_kuuden_henkilon_asuntokunnat_osuus, Vahintaan_seitseman_henkilon_asuntokunnat_osuus)



#Siviilisäädyt 19-49 vuotiaat
nimet_selittajat<-c("Kunta", "Vuosi","Yhteensa","Naimaton","Naimisissa","Eronnut","Leski")
siviilisaadyt<-read.csv2(file="ainestot/2025_aluejako/siviilisaadyt.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_selittajat)

siviilisaadyt <- siviilisaadyt %>%
  #Korjataan Maarianhaminan kirjoitusasu
  mutate(Kunta = if_else(Kunta == "Maarianhamina - Mariehamn", "Maarianhamina", Kunta)
  )

siviilisaadyt <- siviilisaadyt %>%
  mutate(Naimaton_osuus = Naimaton/Yhteensa,
         Naimisissa_osuus = Naimisissa/Yhteensa,
         Eronnut_osuus = Eronnut/Yhteensa,
         Leski_osuus = Leski/Yhteensa
  )%>%
  select(Kunta, Vuosi, Naimaton_osuus, Naimisissa_osuus, Eronnut_osuus, Leski_osuus)

#Koulutusasteet 15-49 vuotiaat
#Korkeakoulutetut
nimet_koulutus<-c("Vuosi", "Kunta","Ika","Koulutusaste", "Korkea_aste")
kk<-read.csv2(file="ainestot/2023_aluejako/korkea_aste_15_49.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_koulutus)

kk<-kk%>%
  select(Vuosi, Kunta, Korkea_aste)

#Vain perusasteen suorittaneet
nimet_koulutus<-c("Vuosi", "Kunta","Ika","Koulutusaste", "Perusaste")
pa<-read.csv2(file="ainestot/2023_aluejako/vain_perusaste_15_49.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_koulutus)

pa<-pa%>%
  select(Vuosi, Kunta, Perusaste)

#Yhteensä
nimet_koulutus<-c("Vuosi", "Kunta","Ika","Koulutusaste", "Yhteensa")
yht<-read.csv2(file="ainestot/2023_aluejako/koulutus_yhteensa_15_49.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_koulutus)

yht<-yht%>%
  select(Vuosi, Kunta, Yhteensa)

kk<-left_join(kk, pa, by = c("Kunta", "Vuosi"))

kk<-left_join(kk, yht, by = c("Kunta", "Vuosi"))

kk <- kk %>%
  #Korjataan Maarianhaminan kirjoitusasu
  mutate(Kunta = if_else(Kunta == "Maarianhamina - Mariehamn", "Maarianhamina", Kunta)
  )

#Tämä aineisto on 2023 aluejaolla, joten muutetaan se 2025 aluejaolle liittämällä Pertunmaa Mäntyharjuun
kuntaliitos_kunnat<-c("Pertunmaa", "Mäntyharju")

pe_ma<-kk%>%
  filter(Kunta %in% kuntaliitos_kunnat)
pe_ma
pe_ma<-pe_ma%>%
  group_by(Vuosi)%>%
  summarise(Kunta = "Mäntyharju", Korkea_aste = sum(Korkea_aste), Perusaste = sum(Perusaste), Yhteensa = sum(Yhteensa))%>%
  ungroup()

kk<-kk%>%
  filter(!Kunta %in% kuntaliitos_kunnat)

kk<-bind_rows(kk, pe_ma)

#Lasketaan suhteelliset osuudet
kk <- kk %>%
  mutate(Perusaste_osuus = Perusaste/Yhteensa,
         Korkea_aste_osuus = Korkea_aste/Yhteensa) %>%
  select(Kunta, Vuosi, Perusaste_osuus, Korkea_aste_osuus)


#Työllisyystilastot 18-64 vuotiaat
nimet_selittajat<-c("Kunta", "Vuosi", "Yhteensa","Tyovoima","Tyolliset","Tyottomat","Tyovoiman_ulkopuolella_olevat","Opiskelijat","Varusmiehet_siviilipalvelusmiehet","Elakelaiset","Muut_tyovoiman_ulkopuolella_olevat")
tt<-read.csv2(file="ainestot/2023_aluejako/tyossakaynti.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_selittajat)

tt <- tt %>%
  #Korjataan Maarianhaminan kirjoitusasu
  mutate(Kunta = if_else(Kunta == "Maarianhamina - Mariehamn", "Maarianhamina", Kunta)
  )

#Tämä aineisto on 2023 aluejaolla, joten muutetaan se 2025 aluejaolle liittämällä Pertunmaa Mäntyharjuun 
pe_ma<-tt%>%
  filter(Kunta %in% kuntaliitos_kunnat)

pe_ma<-pe_ma%>%
  group_by(Vuosi)%>%
  summarise(Kunta = "Mäntyharju", Tyovoima = sum(Tyovoima), Tyolliset = sum(Tyolliset), Yhteensa = sum(Yhteensa), Tyottomat = sum(Tyottomat), Tyovoiman_ulkopuolella_olevat = sum(Tyovoiman_ulkopuolella_olevat), Opiskelijat = sum(Opiskelijat), Varusmiehet_siviilipalvelusmiehet = sum(Varusmiehet_siviilipalvelusmiehet), Elakelaiset = sum(Elakelaiset), Muut_tyovoiman_ulkopuolella_olevat = sum(Muut_tyovoiman_ulkopuolella_olevat))%>%
  ungroup()

tt<-tt%>%
  filter(!Kunta %in% kuntaliitos_kunnat)

tt<-bind_rows(tt, pe_ma)

#Lasketaan suhteelliset osuudet
tt <- tt %>%
  mutate(Tyovoiman_ulkopuolella_osuus = Tyovoiman_ulkopuolella_olevat/Yhteensa,
         Tyottomat_osuus = Tyottomat/Tyovoima,
         Opiskelijat_osuus = Opiskelijat/Yhteensa,
         Muut_tyovoiman_ulkopuolella_olevat_osuus = Muut_tyovoiman_ulkopuolella_olevat/Yhteensa,
         Varusmiehet_siviilipalvelusmiehet_osuus = Varusmiehet_siviilipalvelusmiehet/Yhteensa,
         Elakelaiset_osuus = Elakelaiset/Yhteensa
         )%>%
  select(Kunta, Vuosi, Tyovoiman_ulkopuolella_osuus, Tyottomat_osuus, Opiskelijat_osuus, Muut_tyovoiman_ulkopuolella_olevat_osuus, Varusmiehet_siviilipalvelusmiehet_osuus, Elakelaiset_osuus)


#Kunnan ulkopuolella tyossä käyvät 18-64 vuotiaat
nimet_selittajat<-c("Kunta", "Vuosi", "Yhteensa","Kunnan_ulkopuolella_tyossa")
ut<-read.csv2(file="ainestot/2023_aluejako/tyo_kunnan_ulkopuolella.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_selittajat)

ut <- ut %>%
  #Korjataan Maarianhaminan kirjoitusasu
  mutate(Kunta = if_else(Kunta == "Maarianhamina - Mariehamn", "Maarianhamina", Kunta)
  )

#Tämä aineisto on 2023 aluejaolla, joten muutetaan se 2025 aluejaolle liittämällä Pertunmaa Mäntyharjuun 
pe_ma<-ut%>%
  filter(Kunta %in% kuntaliitos_kunnat)

pe_ma<-pe_ma%>%
  group_by(Vuosi)%>%
  summarise(Kunta = "Mäntyharju", Kunnan_ulkopuolella_tyossa = sum(Kunnan_ulkopuolella_tyossa), Yhteensa = sum(Yhteensa))%>%
  ungroup()

ut<-ut%>%
  filter(!Kunta %in% kuntaliitos_kunnat)

ut<-bind_rows(ut, pe_ma)

#Lasketaan suhteellinen osuus
ut <- ut %>%
  mutate(Kunnan_ulkopuolella_tyossa = Kunnan_ulkopuolella_tyossa/Yhteensa)%>%
  select(Kunta, Vuosi, Kunnan_ulkopuolella_tyossa)


#Muut työvoiman ulkopuolellat olevat sukupuoleittain 18-64 vuotiaat
nimet_selittajat<-c("Kunta", "Vuosi", "Miehet_Yhteensa","Miehet_muut_tyon_ulkopuolella","Naiset_Yhteensa","Naiset_muut_tyon_ulkopuolella")
muut<-read.csv2(file="ainestot/2023_aluejako/muut_tyon_ulkopuolella.csv", fileEncoding = "windows-1252", header=TRUE, col.names = nimet_selittajat)

muut <- muut %>%
  #Korjataan Maarianhaminan kirjoitusasu
  mutate(Kunta = if_else(Kunta == "Maarianhamina - Mariehamn", "Maarianhamina", Kunta)
  )

#Tämä aineisto on 2023 aluejaolla, joten muutetaan se 2025 aluejaolle liittämällä Pertunmaa Mäntyharjuun 
pe_ma<-muut%>%
  filter(Kunta %in% kuntaliitos_kunnat)

pe_ma<-pe_ma%>%
  group_by(Vuosi)%>%
  summarise(Kunta = "Mäntyharju", Miehet_Yhteensa = sum(Miehet_Yhteensa), Miehet_muut_tyon_ulkopuolella = sum(Miehet_muut_tyon_ulkopuolella), Naiset_Yhteensa = sum(Naiset_Yhteensa), Naiset_muut_tyon_ulkopuolella = sum(Naiset_muut_tyon_ulkopuolella))%>%
  ungroup()

muut<-muut%>%
  filter(!Kunta %in% kuntaliitos_kunnat)

muut<-bind_rows(muut, pe_ma)

#Lasketaan suhteelliset osuudet
muut <- muut %>%
  mutate(Miehet_muut_tyon_ulkopuolella_osuus = Miehet_muut_tyon_ulkopuolella/Miehet_Yhteensa,
         Naiset_muut_tyon_ulkopuolella_osuus = Naiset_muut_tyon_ulkopuolella/Naiset_Yhteensa)%>%
  select(Kunta, Vuosi, Miehet_muut_tyon_ulkopuolella_osuus,Naiset_muut_tyon_ulkopuolella_osuus)



#Yhdistetään selittäjät
selittajat<-tunnuslukuja
selittajat<-left_join(selittajat, kk , by = c("Kunta", "Vuosi")) 
selittajat<-left_join(selittajat, tt, by = c("Kunta", "Vuosi")) 
selittajat<-left_join(selittajat, ut, by = c("Kunta", "Vuosi"))
selittajat<-left_join(selittajat, muut, by = c("Kunta", "Vuosi"))
selittajat<-left_join(selittajat, siviilisaadyt, by = c("Kunta", "Vuosi")) 
selittajat<-left_join(selittajat, asuntokunnat, by = c("Kunta", "Vuosi"))
selittajat<-left_join(selittajat, naiset_19_49, by = c("Kunta", "Vuosi")) 

#Kasvatetaan selittäjissä vuotta yhdellä, koska tiedot ovat vuoden viimeiseltä päivältä
selittajat <- selittajat %>%
  mutate(Vuosi = Vuosi + 1)

#Puuttuvat tiedot
selittajat[selittajat == "."] <- NA

#Pudotetaan vuosi 2025, koska sitä ei ole vasteessa
selittajat<-selittajat%>%
  filter(Vuosi<2025)

#Skaalataan selittajat ja muutetaan liukuluvuiksi
selittajat_skaalattu <- selittajat %>%
  mutate(across(!all_of(c("Kunta", "Vuosi")), ~ scale(as.double(.))))

selittajat <- selittajat %>%
  mutate(across(!all_of(c("Kunta", "Vuosi")), ~ as.double(.)))
 
#Liitetään selittäjät päätauluun
kunnat_skaalattu<-left_join(kunnat, selittajat_skaalattu, by = c("Kunta", "Vuosi")) 
kunnat<-left_join(kunnat, selittajat, by = c("Kunta", "Vuosi"))

#Vuosittaiset kokonaishedelmällisyyslukujen keskiarvot
keskiarvot_hed_vuosittain <- kunnat %>%
  group_by(Vuosi) %>%
  summarize(Hedelmallisyys = mean(Hedelmallisyys))

head(kunnat_skaalattu)
#-------------------------------------------
#Puuttuva tieto
#-------------------------------------------
sum(is.na(kunnat%>%filter(Vuosi>1990)))

#-------------------------------------------
# Tunnusluvut
#-------------------------------------------
#Suurin kokonaishedelmällisyysluku
kunnat[which.max(kunnat$Hedelmallisyys),c("Vuosi", "Kunta", "Hedelmallisyys")]

#Havaintojen lukumäärä, joilla kokonaishedelmällisyysluku on nolla
dim(kunnat[which(kunnat$Hedelmallisyys==0),c("Vuosi", "Kunta", "Hedelmallisyys")])

#Kokonaishedelmällisyyslukujen keskiarvo ja keskihajonta
mean(kunnat$Hedelmallisyys)
sd(kunnat$Hedelmallisyys)

#Havainnot, joilla poikkeavan suuret arvot
kunnat %>%
  filter(Hedelmallisyys > 6) %>%
    select(Kunta, Vuosi, Hedelmallisyys)

#-------------------------------------------
# Visualisointi
#-------------------------------------------
#Suomen kokonaishedelmällisyys
suomen_hedelmallisyys<-read.csv2(file="ainestot/suomen_hedelmallisyysluku.csv", fileEncoding = "windows-1252", header=TRUE, col.names = c("Vuosi","Kokonaishedelmällisyysluku"))

suomen_hedelmallisyys <- suomen_hedelmallisyys %>%
  mutate(Kokonaishedelmällisyysluku = as.double(Kokonaishedelmällisyysluku))

suomen_hedelmallisyys_plot <- ggplot(suomen_hedelmallisyys, aes(y=Kokonaishedelmällisyysluku, x=Vuosi)) + 
  geom_line() +
  ylim(
    c(0,6)
  )
suomen_hedelmallisyys_plot

#ggsave("suomen_hedelmallisyys.pdf",suomen_hedelmallisyys_plot, width=16, height=10, units="cm")

#Kokonaishedelmällisyyslukujen histogrammi
kunnat_hist <- ggplot(kunnat, aes(x=Hedelmallisyys)) + 
  geom_histogram(binwidth = 0.05, alpha = 1) +
  labs(y = "Lukumäärä", x = "Kokonaishedelmällisyysluku")

kunnat_hist
#ggsave("kunnat_hist.pdf",kunnat_hist, width=9, height=9, units="cm")


#Kokonaishedelmällisyyslukujen keskiarvo vuosittain
kunnat_aika <- ggplot(keskiarvot_hed_vuosittain, aes(y=Hedelmallisyys, x=Vuosi)) + 
  geom_line() +
  labs(y = "Kokonaishedelmällisyysluku", x = "Vuosi")

kunnat_aika
#ggsave("kunnat_aika.pdf",kunnat_aika, width=9, height=9, units="cm")

#Suurin vuosittainen keskiarvo
keskiarvot_hed_vuosittain[which.max(keskiarvot_hed_vuosittain$Hedelmallisyys),]

#Pienin vuosittainen keskiarvo
keskiarvot_hed_vuosittain[which.min(keskiarvot_hed_vuosittain$Hedelmallisyys),]

#-------------------
# Paikkatiedot
#-------------------
#Kokonaishedelmällisyys alueellisesti
polygon1 <- get_municipalities(year = 2025, scale = 4500)

#Liitetään maatieteelliset tiedot päätauluun
polygon_join <- polygon1 %>%
  select(kunta, maakunta_name_fi)

kunnat_skaalattu<-left_join(kunnat_skaalattu, polygon_join, by=join_by(Id == kunta))

kunnat_skaalattu<-kunnat_skaalattu%>%
  mutate(Maakunta = as.factor(maakunta_name_fi))%>%
  select(!geom)

#Lasketaan keskiarvot
keskiarvot <- kunnat %>%
  group_by(Id) %>%
  summarise(Hedelmallisyys = mean(Hedelmallisyys))

keskiarvot_p<-right_join(polygon1, keskiarvot, by=join_by(kunta == Id))

#Suurin keskiarvo
keskiarvot_p[which.max(keskiarvot_p$Hedelmallisyys),c("Hedelmallisyys", "nimi")]

#Pienin keskiarvo
keskiarvot_p[which.min(keskiarvot_p$Hedelmallisyys),c("Hedelmallisyys", "nimi")]

#Karttakuva kuntien keskiarvoista vuosilta 1990 - 2024
kunnat_alue <- ggplot(keskiarvot_p, 
  aes(fill = Hedelmallisyys)) +
  geom_sf() +
  scale_fill_gradient(low = "white", high = "red", name ="Kokonais-\nhedelmällisyysluku") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

kunnat_alue
#ggsave("kunnat_alue.pdf",kunnat_alue, width=12, height=12, units="cm")

#Maakuntien kuntien kokonaishedelmällisyyslukujen ajallinen vaihtelu
maakunnat_keskiarvot <- kunnat_skaalattu %>%
  group_by(Maakunta, Vuosi) %>%
  summarise(ka_Hedelmallisyys = mean(Hedelmallisyys))

maakunnat_koko_ka <- maakunnat_keskiarvot %>%
  group_by(Maakunta) %>%
  summarise(koko_ka_Hedelmallisyys = mean(ka_Hedelmallisyys))

maakunnat_yht <- left_join(maakunnat_keskiarvot, maakunnat_koko_ka, by = "Maakunta")

maakunnat_yht <- maakunnat_yht %>%
  mutate(ero = abs(ka_Hedelmallisyys-koko_ka_Hedelmallisyys)) %>%
  group_by(Maakunta) %>%
  summarize(max_ero = max(ero))

maakunnat_yht <- maakunnat_yht %>%
  arrange(max_ero)

maakunnat_yht

#--------------------------------------
# Mallinnus
#--------------------------------------
#Korvataan aineistosta 0 arvolla 0.1
ei_nollia_kunnat_skaalattu<-kunnat_skaalattu%>%
  mutate(Hedelmallisyys = if_else(Hedelmallisyys==0, 0.1, Hedelmallisyys))

#Asetetaan refrenssiluokaksi Pohjanmaa
ei_nollia_kunnat_skaalattu$Maakunta <- relevel(ei_nollia_kunnat_skaalattu$Maakunta, ref = "Pohjanmaa")

#Poistetaan turhat muuttujat
rm(naiset_19_49, selittajat, selittajat_skaalattu, syntyneet, tunnuslukuja, vakiluku, yht, pa, kk, tt, ut, muut, siviilisaadyt, asuntokunnat, kunnat_afr, pe_ma)

#--------------
# Normaalijakaumamallit
#--------------
#
#Yhteinen varianssi
#
malli_norm<-obs(Hedelmallisyys ~ -1 + varying(~ 1 + lag(Hedelmallisyys, k=1) + Maakunta + Korkea_aste_osuus + Perusaste_osuus + Tyottomat_osuus + Miehet_muut_tyon_ulkopuolella_osuus + Opiskelijat_osuus + Naiset_muut_tyon_ulkopuolella_osuus + Kunnan_ulkopuolella_tyossa + log_Vakiluku + Ei_uskontokuntaa + Keski_ika + Asuinalueella_syntyneet + Naimisissa_osuus + Yhden_henkilon_asuntokunnat_osuus + Neljan_henkilon_asuntokunnat_osuus + Vahintaan_seitseman_henkilon_asuntokunnat_osuus + Naiset_osuus_19_49) + random(~ 1),
                family = "gaussian") + splines(df = 10, lb_tau = 0.01)

fit_malli_norm <- dynamite(
  dformula = malli_norm,
  data = ei_nollia_kunnat_skaalattu, time = "Vuosi", group = "Kunta",
  chains = 4, cores = 4, seed = 189, refresh = 1, warmup = 10000, iter = 30000, thin = 10
)

#Tallennetaan malli
#saveRDS(fit_malli_norm, "fit_malli_norm_lag.rds")
fit_malli_norm <- readRDS("fit_malli_norm_lag.rds")

#Konvergenssin tarkastelu
ess_r(fit_malli_norm)

#Sovitteen posteriorinäytteet
sovitteet_norm <- fitted(fit_malli_norm, thin = 4)

#Tallennetaan sovitteet
saveRDS(sovitteet_norm, "sovitteet_norm.rds")
#sovitteet_norm <- readRDS("sovitteet_norm.rds")

#R2
norm_r2 <-R2(sovitteet_norm)

#Posteriorikeskiarvo
mean(norm_r2)

#95% posterioriväli
quantile(norm_r2, probs = c(0.025, 0.975))

#
#Kuntakohtainen varianssi, joka mallinnettu regressiomallilla
#
#Stan-koodi
uusi_stan_koodi_norm<-"functions {
}
data {
  int<lower=1> T; // number of time points
  int<lower=1> N; // number of individuals
  int<lower=0> K; // total number of covariates across all channels
  array[T] matrix[N, K] X; // covariates as an array of N x K matrices
  row_vector[K] X_m; // Means of all covariates at first time point
  int<lower=1> D; // number of B-splines
  matrix[D, T] Bs; // B-spline basis matrix
  int<lower=0> M; // number group-level effects (including intercepts)
  // number of fixed, varying and random coefficients, and related indices
  int<lower=0> K_varying_Hedelmallisyys;
  int<lower=0> K_random_Hedelmallisyys; // includes the potential random intercept
  int<lower=0> K_Hedelmallisyys; // K_fixed + K_varying
  array[K_varying_Hedelmallisyys] int J_varying_Hedelmallisyys;
  array[K_Hedelmallisyys] int J_Hedelmallisyys; // fixed and varying
  array[K_varying_Hedelmallisyys] int L_varying_Hedelmallisyys;
  // Parameters of vectorized priors
  matrix[K_varying_Hedelmallisyys, 2] delta_prior_pars_Hedelmallisyys;
  matrix[K_varying_Hedelmallisyys, 2] tau_prior_pars_Hedelmallisyys;
  matrix[K_random_Hedelmallisyys, 2] sigma_nu_prior_pars_Hedelmallisyys;
  matrix[N, T] y_Hedelmallisyys;
}
transformed data {
  matrix[N,T] log_V;
  matrix[N,1] log_V_mean;
  for (t in 1:T){
    log_V[,t] = X[t][, 27];
  }
  for (n in 1:N){
    log_V_mean[n,1] = mean(log_V[n,]);
  }
}
parameters {
  // Random group-level effects
  vector<lower=0>[M] sigma_nu; // standard deviations of random effects
  matrix[N, M] nu_raw;
  matrix[K_varying_Hedelmallisyys, D] omega_Hedelmallisyys; // Spline coefficients
  vector<lower=0.01>[K_varying_Hedelmallisyys] tau_Hedelmallisyys; // SDs for the random walks
  real a_Hedelmallisyys; // Mean of the first time point
  row_vector[D - 1] omega_raw_alpha_Hedelmallisyys; // Coefficients for alpha
  real<lower=0.01> tau_alpha_Hedelmallisyys; // SD for the random walk
  vector<lower=0>[N] sigma_Hedelmallisyys; // SD of the normal distribution
  real<lower=0> sigma_sigma; // SD of the sigma prior
  real alpha_sigma; // Alpha parameter of the sigma prior
  vector[1] beta_sigma; // Beta parameter of the sigma prior
}
transformed parameters {
  vector[1] sigma_nu_Hedelmallisyys = sigma_nu[1:1];
  matrix[N, 1] nu_Hedelmallisyys = diag_post_multiply(nu_raw[, 1:1], sigma_nu_Hedelmallisyys);
  // Time-varying coefficients
  array[T] vector[K_varying_Hedelmallisyys] delta_Hedelmallisyys;
  // Time-varying intercept
  array[T] real alpha_Hedelmallisyys;
  // Spline coefficients
  real omega_alpha_1_Hedelmallisyys;
  row_vector[D] omega_alpha_Hedelmallisyys;
  for (t in 1:T) {
    delta_Hedelmallisyys[t] = omega_Hedelmallisyys * Bs[, t];
  }
  // Define the first alpha using mean a_Hedelmallisyys
  {
    vector[K_Hedelmallisyys] gamma__Hedelmallisyys;
    gamma__Hedelmallisyys[L_varying_Hedelmallisyys] = delta_Hedelmallisyys[1];
    omega_alpha_1_Hedelmallisyys = a_Hedelmallisyys - X_m[J_Hedelmallisyys] * gamma__Hedelmallisyys;
  }
  omega_alpha_Hedelmallisyys[1] = omega_alpha_1_Hedelmallisyys;
  omega_alpha_Hedelmallisyys[2:D] = omega_raw_alpha_Hedelmallisyys;
  for (t in 1:T) {
    alpha_Hedelmallisyys[t] = omega_alpha_Hedelmallisyys * Bs[, t];
  }
}
model {
  to_vector(nu_raw) ~ std_normal();
  sigma_nu_Hedelmallisyys ~ normal(sigma_nu_prior_pars_Hedelmallisyys[, 1], sigma_nu_prior_pars_Hedelmallisyys[, 2]);
  a_Hedelmallisyys ~ normal(2, 2);
  omega_raw_alpha_Hedelmallisyys[1] ~ normal(omega_alpha_1_Hedelmallisyys, tau_alpha_Hedelmallisyys);
  for (i in 2:(D - 1)) {
    omega_raw_alpha_Hedelmallisyys[i] ~ normal(omega_raw_alpha_Hedelmallisyys[i - 1], tau_alpha_Hedelmallisyys);
  }
  tau_alpha_Hedelmallisyys ~ normal(0, 2);
  omega_Hedelmallisyys[, 1] ~ normal(delta_prior_pars_Hedelmallisyys[, 1], delta_prior_pars_Hedelmallisyys[, 2]);
  for (i in 2:D) {
    omega_Hedelmallisyys[, i] ~ normal(omega_Hedelmallisyys[, i - 1], tau_Hedelmallisyys);
  }
  tau_Hedelmallisyys ~ normal(tau_prior_pars_Hedelmallisyys[, 1], tau_prior_pars_Hedelmallisyys[, 2]);
  sigma_sigma ~ exponential(1);
  alpha_sigma ~ normal(0, 2);
  beta_sigma ~ normal(0, 2);
  log(sigma_Hedelmallisyys) ~ normal_id_glm(log_V_mean, alpha_sigma, beta_sigma, sigma_sigma);
  {
    real ll = 0.0;
    vector[K_Hedelmallisyys] gamma__Hedelmallisyys;
    for (t in 1:T) {
      vector[N] intercept_Hedelmallisyys = alpha_Hedelmallisyys[t] + nu_Hedelmallisyys[, 1];
      gamma__Hedelmallisyys[L_varying_Hedelmallisyys] = delta_Hedelmallisyys[t];
      ll += normal_id_glm_lupdf(y_Hedelmallisyys[, t] | X[t][, J_Hedelmallisyys], intercept_Hedelmallisyys, gamma__Hedelmallisyys, sigma_Hedelmallisyys);
    }
    target += ll;
  }
}
generated quantities {
}
"

fit_malli_norm_var <- dynamite(
  dformula = malli_norm,
  data = ei_nollia_kunnat_skaalattu, 
  custom_stan_model = uusi_stan_koodi_norm,
  time = "Vuosi", group = "Kunta",
  chains = 4, cores = 4, seed = 189, refresh = 1, iter = 30000, warmup = 10000, thin=10
)

#Tallennetaan malli
#saveRDS(fit_malli_norm_var, "fit_malli_norm_var_lag.rds")
fit_malli_norm_var <- readRDS("fit_malli_norm_var_lag.rds")

#Konvergenssin tarkastelu
ess_r(fit_malli_norm_var)

#Parametrien tulkinta

#Alphan posteriorikeskiarvot ja 95% posteriorivälit
alpha <- coef(fit_malli_norm_var, types = "alpha", probs = c(0.025, 0.975))
print(alpha, n = 35)

#Posteriorinäytteet maakunnan delta-parametreille
naytteet_delta_maakunnat <- as_draws(fit_malli_norm_var, parameters = get_parameter_names(fit_malli_norm_var)[3:20])

#Lasketaan tunnusluvut parametreille
delta_maakunnat_tunnusluvut <- posteriori_tunnusluvut(naytteet_delta_maakunnat)

#Nollasta merkittävästi poikkeavat maakunnat
#Etelä-Karjala
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Etelä-Karjala")) %>%
  print(n=35)

#Etelä-Pohjanmaa
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Etelä-Pohjanmaa")) %>%
  print(n=35)

#Etelä-Savo
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Etelä-Savo")) %>%
  print(n=35)

#Kainuu
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Kainuu")) %>%
  print(n=35)

#Kanta-Häme
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Kanta-Häme")) %>%
  print(n=35)

#Keski-Pohjanmaa
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Keski-Pohjanmaa")) %>%
  print(n=35)

#Keski-Suomi
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Keski-Suomi")) %>%
  print(n=35)

#Kymenlaakso
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Kymenlaakso")) %>%
  print(n=35)

#Lappi
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Lappi")) %>%
  print(n=35)

#Pirkanmaa
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Pirkanmaa")) %>%
  print(n=35)

#Pohjois-Karjala
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Pohjois-Karjala")) %>%
  print(n=35)

#Pohjois-Pohjanmaa
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Pohjois-Pohjanmaa")) %>%
  print(n=35)

#Pohjois-Savo
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Pohjois-Savo")) %>%
  print(n=35)

#Päijät-Häme
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Päijät-Häme")) %>%
  print(n=35)

#Satakunta
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Satakunta")) %>%
  print(n=35)

#Uusimaa
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Uusimaa")) %>%
  print(n=35)

#Varsinais-Suomi
delta_maakunnat_tunnusluvut %>%
  filter(str_detect(Selittaja, "Varsinais-Suomi")) %>%
  print(n=35)

#Posteriorinäytteet muille delta-parametreille
naytteet_delta_muut <- as_draws(fit_malli_norm_var, parameters = get_parameter_names(fit_malli_norm_var)[c(2,21:36)])

delta_muut_tunnusluvut <- posteriori_tunnusluvut(naytteet_delta_muut)

#Asuinalueellaan syntyneiden osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Asuinalueella")) %>%
  print(n=35)

#Uskontokuntaan kuulumattomien osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Ei_uskontokuntaa")) %>%
  print(n=35)

#Edellisen vuoden kokonaishedelmällisyysluku
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "lag1")) %>%
  print(n=35)

#Keski-ikä
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Keski_ika")) %>%
  print(n=35)

#Korkeakoulutettujen osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Korkea")) %>%
  print(n=35)

#Kunnan ulkopuolella työssä käyvien osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Kunnan_ulkopuolella")) %>%
  print(n=35)

#Väkiluvun logaritmi
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Vakiluku")) %>%
  print(n=35)

#Muiden työn ulkopuolella olevien osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Miehet_muut")) %>%
  print(n=35)

#Naimisissa olevien osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Naimisissa")) %>%
  print(n=35)

#Muiden työn ulkopuolella olevien osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Naiset_muut")) %>%
  print(n=35)

#Naisten osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Naiset_osuus")) %>%
  print(n=35)

#Neljän henkilön asuntokuntien osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Neljan")) %>%
  print(n=35)

#Opiskelijoiden osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Opiskeli")) %>%
  print(n=35)

#Vain perusasteen suorittaneiden osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Perusaste")) %>%
  print(n=35)

#Työttömien osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Tyotto")) %>%
  print(n=35)

#Vähintään seitsemän henkilön asuntokuntien osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "seitseman")) %>%
  print(n=35)

#Yhden henkilön asuntokuntien osuus
delta_muut_tunnusluvut %>%
  filter(str_detect(Selittaja, "Yhden")) %>%
  print(n=35)


#Kuvien otsikot
facet_alpha <- c("alpha")
names(facet_alpha) <- c("alpha_Hedelmallisyys")

maakunta_selittaja_nimet <-  get_parameter_names(fit_malli_norm_var)[3:20]

facet_delta_maakunnat <- str_remove(maakunta_selittaja_nimet, "delta_Hedelmallisyys_Maakunta")
names(facet_delta_maakunnat) <- maakunta_selittaja_nimet

facet_delta_muut <- c("Edellisen vuoden \n kokonaishedelmällisyysluku", "Korkeakoulutettujen osuus", "Vain perusasteen suorittaneiden \n osuus", "Työttömien osuus"                                
                      , "Muiden työn ulkopuolella olevien \n miesten osuus", "Opiskelijoiden osuus", "Muiden työn ulkopuolella olevien \n naisten osuus", "Kunnan ulkopuolella \n työssäkäyvien osuus"
                      , "log(Väkiluku)", "Uskontokuntaan kuulumattomien \n osuus", "Keski-ikä", "Asuinalueellaan syntyneiden \n osuus"
                      , "Naimisissa olevien osuus", "Yhden henkilön asuntokuntien \n osuus", "Neljän henkilön asuntokuntien \n osuus", "Vähintään seitsemän henkilön \n asuntokuntien osuus"
                      , "Naisten osuus")
names(facet_delta_muut) <- get_parameter_names(fit_malli_norm_var)[c(2,21:36)]

#Kuvat selittäjien posteriorikeskiarvoista ja 95% posterioriväleistä
#Alpha
p_alpha <- plot(fit_malli_norm_var, type=c("alpha"), n_params = 1000, level = 0.025, scales = "fixed") + 
  ggtitle(element_blank()) + xlab("Vuosi") + ylab("Arvo") + 
  facet_wrap("parameter",labeller = labeller(parameter = facet_alpha))
p_alpha
ggsave("alpha.pdf", p_alpha, width=8, height=8, units="cm")

#Maakunta deltat
p_delta_maakunnat <- plot(fit_malli_norm_var, type=c("delta"), n_params = 1000, parameters = names(facet_delta_maakunnat), level = 0.025) + 
  ggtitle(element_blank()) + xlab("Vuosi") + ylab("Arvo") + 
  facet_wrap("parameter",labeller = labeller(parameter = facet_delta_maakunnat), ncol = 4, scales = "fixed") + geom_hline(yintercept = 0, linetype="dashed", color = "red")
p_delta_maakunnat
ggsave("delta_maakunnat.pdf",p_delta_maakunnat, width=16, height=18, units="cm")

#Muut deltat
p_delta_muut <- plot(fit_malli_norm_var, type=c("delta"), n_params = 1000, parameters = names(facet_delta_muut), level = 0.025) + 
  ggtitle(element_blank()) + xlab("Vuosi") + ylab("Arvo") + 
  facet_wrap("parameter",labeller = labeller(parameter = facet_delta_muut), ncol = 3, scales = "free") + geom_hline(yintercept = 0, linetype="dashed", color = "red")
p_delta_muut
ggsave("delta_muut.pdf",p_delta_muut, width=19, height=25, units="cm")

#Esitykseen kuva
#Muut deltat
p_delta_muut_esitys <- plot(fit_malli_norm_var, type=c("delta"), n_params = 1000, parameters = names(facet_delta_muut), level = 0.025) + 
  ggtitle(element_blank()) + xlab("Vuosi") + ylab("Arvo") + 
  facet_wrap("parameter",labeller = labeller(parameter = facet_delta_muut), ncol = 4, scales = "free") + geom_hline(yintercept = 0, linetype="dashed", color = "red")
p_delta_muut_esitys
#ggsave("delta_muut_esitys.pdf",p_delta_muut_esitys, width=25, height=19, units="cm")

#Sigma_nu:n posteriorijakauma
sigma_nu <- coef(fit_malli_norm_var, types = "sigma_nu", probs = c(0.025, 0.975))
sigma_nu

#Satunnaisvaikutusten visualisointi
nu <- coef(fit_malli_norm_var, types = "nu", probs = c(0.025, 0.975))

nu$group <- factor(nu$group, levels = nu$group[order(nu$mean)])

nu_poikkeamat <- nu %>%
  filter(q2.5 > 0 | q97.5 < 0)

nu_poikkeamat <- nu_poikkeamat %>%
  arrange(mean)

nu_poikkeamat
p_nu <- ggplot(data = nu_poikkeamat, aes(x = mean, y = group)) +
  geom_pointrange(aes(x = mean, y = group, xmin = q2.5, xmax = q97.5), fatten = 1) +
  geom_vline(xintercept = 0, linetype="dashed", color = "red") +
  labs(x = "Arvo", y = "Kunta")
p_nu
#ggsave("nu_poikkeamat.pdf",p_nu, width=12, height=16, units="cm")

#Aineisto karttakuvaa varten
kunnat_nu <- kunnat_skaalattu %>%
  filter(Vuosi == 2000) %>%
  select(Kunta, Id)

kunnat_nu <- left_join(kunnat_nu, nu[,c("group", "mean")], by = join_by(Kunta == group))

nu_kartta <- right_join(polygon1, kunnat_nu, by=join_by(kunta == Id))

#Karttakuva satunnaisvaikutusten posteriorikeskiarvoille
nu_kuva <- ggplot(nu_kartta, aes(fill = mean)) +
  geom_sf() +
  scale_fill_gradient(low = "white", high = "red", name ="Satunnaisvaikutusten\nposteriorikeskiarvo") +
  theme(axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_blank())

nu_kuva
#ggsave("nu_kartta.pdf",nu_kuva, width=12, height=12, units="cm")

#Varianssin posteriorinäytteet ja niiden keskiarvot
naytteet<-posterior::as_draws_df(fit_malli_norm_var$stanfit)

#Sigma^(sigma), alpha^(sigma) ja beta^(sigma)-parametrien posteriorinäytteet
naytteet_sigma <- data.frame(Arvo = naytteet$sigma_sigma, Muuttuja = "sigma_sigma")
naytteet_sigma <- rbind(naytteet_sigma, data.frame(Arvo = naytteet$alpha_sigma, Muuttuja = "alpha_sigma"))
naytteet_sigma <- rbind(naytteet_sigma, data.frame(Arvo = naytteet$`beta_sigma[1]`, Muuttuja = "beta_sigma"))

#Tunnusluvut
naytteet_sigma_kpi <- naytteet_sigma%>% 
  group_by(Muuttuja)%>%
  summarize(Keskiarvo = round(mean(Arvo),4), Ala95 = quantile(Arvo, probs = 0.025), Yla95 = quantile(Arvo, probs = 0.975)) %>%
  mutate(across(all_of(c("Keskiarvo", "Ala95", "Yla95")), ~ format(round(.,digits=3), nsmall = 3)))
  
naytteet_sigma_kpi

#Muutetaan dynamiten prepare_eval_env_univariate-funktiota, jotta predict-funktio
#toimii yksilökohtaisilla hajonta- ja muotoparametreilla
prepare_eval_env_univariate_oma <- function (e, resp, resp_levels, cvars, samples, nu_samples, has_random_effects, 
                                             idx, type, eval_type) 
{
  alpha <- paste0("alpha_", resp)
  beta <- paste0("beta_", resp)
  cutpoint <- paste0("cutpoint_", resp)
  delta <- paste0("delta_", resp)
  phi <- paste0("phi_", resp)
  sigma <- paste0("sigma_", resp)
  lambda <- paste0("lambda_", resp)
  psi <- paste0("psi_", resp)
  e$J_fixed <- cvars$J_fixed
  e$K_fixed <- cvars$K_fixed
  e$J_varying <- cvars$J_varying
  e$K_varying <- cvars$K_varying
  e$J_random <- cvars$J_random
  e$K_random <- cvars$K_random
  e$has_random_intercept <- cvars$has_random_intercept
  e$has_lfactor <- cvars$has_lfactor
  e$resp <- resp
  e$phi <- onlyif(!is.null(samples[[phi]]), 
                    if(length(c(samples[[phi]][idx, , drop = FALSE])) == e$k){
                      c(samples[[phi]][idx, , drop = FALSE])
                    }
                    else{
                      rep_len(c(samples[[phi]][idx]), length.out = e$k)
                    })
  e$sigma <- onlyif(!is.null(samples[[sigma]]), 
                    if(length(c(samples[[sigma]][idx, , drop = FALSE])) == e$k){
                      c(samples[[sigma]][idx, , drop = FALSE])
                    }
                    else{
                      rep_len(c(samples[[sigma]][idx]), length.out = e$k)
                    })

  if (has_random_effects) {
    nus <- make.unique(rep(paste0("nu_", resp), e$K_random))
    e$nu <- nu_samples[, , nus, drop = FALSE]
  }
  if (is_cumulative(e$family)) {
    e$d <- cvars$S
    e$mean_cols <- paste0(resp, "_mean_", resp_levels)
    e$fitted_cols <- paste0(resp, "_fitted_", resp_levels)
    e$invlink <- ifelse_(identical(e$family$link, "logit"), 
                         stats::plogis, stats::pnorm)
    if (cvars$has_fixed_intercept) {
      e$cutpoint <- samples[[cutpoint]][idx, , drop = FALSE]
      e$cutpoint <- e$cutpoint[rep_len(e$n_draws, e$k), 
                               , drop = FALSE]
      e$alpha <- matrix(0, e$n_draws, 1L)
    }
    if (cvars$has_varying_intercept) {
      e$cutpoint <- samples[[cutpoint]][idx, , , drop = FALSE]
      e$cutpoint <- e$cutpoint[rep_len(e$n_draws, e$k), 
                               , , drop = FALSE]
      e$alpha <- matrix(0, e$n_draws, dim(e$cutpoint)[2L])
    }
  }
  else {
    if (cvars$has_fixed_intercept) {
      e$alpha <- array(samples[[alpha]][idx], c(e$n_draws, 
                                                1L))
    }
    if (cvars$has_varying_intercept) {
      e$alpha <- samples[[alpha]][idx, , drop = FALSE]
    }
  }
  if (cvars$has_lfactor) {
    e$lambda <- samples[[lambda]][idx, , drop = FALSE]
    e$psi <- samples[[psi]][idx, , drop = FALSE]
  }
  e$beta <- samples[[beta]][idx, , drop = FALSE]
  e$delta <- samples[[delta]][idx, , , drop = FALSE]
  e$xbeta <- numeric(e$k)
  e$call <- generate_sim_call_univariate(resp = resp, family = e$family, 
                                         type = type, eval_type = eval_type, has_fixed = cvars$has_fixed, 
                                         has_varying = cvars$has_varying, has_random = cvars$has_random, 
                                         has_fixed_intercept = cvars$has_fixed_intercept, has_varying_intercept = cvars$has_varying_intercept, 
                                         has_random_intercept = cvars$has_random_intercept, has_offset = cvars$has_offset, 
                                         has_lfactor = cvars$has_lfactor)
}

onlyif <- get("onlyif", envir = asNamespace("dynamite"))
is_cumulative <- get("is_cumulative", envir = asNamespace("dynamite"))
generate_sim_call_univariate <- get("generate_sim_call_univariate", envir = asNamespace("dynamite"))

assignInNamespace("prepare_eval_env_univariate", prepare_eval_env_univariate_oma, ns = "dynamite")

#Ennustejakauman posteriorinäytteet
set.seed(189)
ennustukset_norm_var <- predict(fit_malli_norm_var, expand = TRUE, thin = 4)

#Sovitteen posteriorinäytteet
sovitteet_norm_var <- fitted(fit_malli_norm_var, thin = 4)

#Tallennetaan ennustukset
#saveRDS(ennustukset_norm_var, "ennustukset_norm_var.rds")
ennustukset_norm_var <- readRDS("ennustukset_norm_var.rds")

#Tallennetaan sovitteet
saveRDS(sovitteet_norm_var, "sovitteet_norm_var.rds")
#sovitteet_norm_var <- readRDS("sovitteet_norm_var.rds")

#Ennusteposteriorinäytteiden ero oikeaan arvoon
ennustukset_norm_var <- ennustukset_norm_var %>%
  filter(Vuosi>1990) %>%
    mutate(Hedelmallisyys_ero = abs(Hedelmallisyys-Hedelmallisyys_new))

#Kaikki erot summattuna
ennustukset_norm_var_erot <- ennustukset_norm_var %>%
  group_by(Kunta) %>%
  summarise(Hedelmallisyys_ero_summa = sum(Hedelmallisyys_ero))

#Keskimääräiset erot per näyte
ennustukset_norm_var_erot <- ennustukset_norm_var_erot %>%
  mutate(Hedelmallisyys_ero_ka = Hedelmallisyys_ero_summa/(length(unique(ennustukset_norm_var$Vuosi))*length(unique(ennustukset_norm_var$.draw))))

ennustukset_norm_var_erot <- ennustukset_norm_var_erot %>%
  arrange(Hedelmallisyys_ero_summa)

#9 parasta ja 9 huonointa
ennustukset_norm_var_erot[1:9,]
ennustukset_norm_var_erot[299:307,]

#Kuvat kunnista, joilla pienin ero ennusteiden ja oikean arvon välillä
ennusteet_tarkat <- ggplot(data = ennustukset_norm_var %>% filter(Kunta %in% as.vector(ennustukset_norm_var_erot$Kunta[1:9])), aes(x = Vuosi)) + 
  stat_summary(aes(y = Hedelmallisyys_new), fun.min = function(y) quantile(y,probs=c(0.025)), fun.max = function(y) quantile(y,probs=c(0.975)), geom="ribbon" , fill = "black", alpha = 0.1) +
  stat_summary(aes(y = Hedelmallisyys_new), fun.min = function(y) quantile(y,probs=c(0.125)), fun.max = function(y) quantile(y,probs=c(0.875)), geom="ribbon" , fill = "black", alpha = 0.2) +
  stat_summary(aes(y = Hedelmallisyys_new), fun.min = function(y) quantile(y,probs=c(0.25)), fun.max = function(y) quantile(y,probs=c(0.75)), geom="ribbon" , fill = "black", alpha = 0.3) +
  geom_line(aes(y = Hedelmallisyys), colour="red", linewidth = 0.7) +
  ylab("Kokonaishedelmällisyysluku") +
  stat_summary(aes(y = Hedelmallisyys_new), fun="mean", colour = "black", geom="line", linewidth = 0.7) +
  facet_wrap(~Kunta, ncol = 3, scales = "free")

ennusteet_tarkat
ggsave("ennusteet_tarkat.pdf",ennusteet_tarkat, width=16, height=16, units="cm")

ennusteet_epatarkat <- ggplot(data = ennustukset_norm_var %>% filter(Kunta %in% as.vector(ennustukset_norm_var_erot$Kunta[299:307])), aes(x = Vuosi)) + 
  stat_summary(aes(y = Hedelmallisyys_new), fun.min = function(y) quantile(y,probs=c(0.025)), fun.max = function(y) quantile(y,probs=c(0.975)), geom="ribbon" , fill = "black", alpha = 0.1) +
  stat_summary(aes(y = Hedelmallisyys_new), fun.min = function(y) quantile(y,probs=c(0.125)), fun.max = function(y) quantile(y,probs=c(0.875)), geom="ribbon" , fill = "black", alpha = 0.2) +
  stat_summary(aes(y = Hedelmallisyys_new), fun.min = function(y) quantile(y,probs=c(0.25)), fun.max = function(y) quantile(y,probs=c(0.75)), geom="ribbon" , fill = "black", alpha = 0.3) +
  geom_line(aes(y = Hedelmallisyys), colour="red", linewidth = 0.7) +
  ylab("Kokonaishedelmällisyysluku") +
  stat_summary(aes(y = Hedelmallisyys_new), fun="mean", colour = "black", geom="line", linewidth = 0.7) +
  facet_wrap(~Kunta, ncol = 3, scales = "free")

ennusteet_epatarkat
ggsave("ennusteet_epatarkat.pdf",ennusteet_epatarkat, width=16, height=16, units="cm")

#Negatiivisten ennusteiden osuus
ennustukset_norm_var %>%
  mutate(Hedelmallisyys_new_nega = ifelse(Hedelmallisyys_new<0, 1, 0)) %>%
  summarise(Negatiiviset_osuus = mean(Hedelmallisyys_new_nega))

#R2
norm_var_r2 <-R2(sovitteet_norm_var)

#Posteriorikeskiarvo
mean(norm_var_r2)

#95% posterioriväli
quantile(norm_var_r2, probs = c(0.025, 0.975))

#
#Kuntakohtainen varianssi itsenäisenä parametrina
#

uusi_stan_koodi_norm_its<-"functions {
}
data {
  int<lower=1> T; // number of time points
  int<lower=1> N; // number of individuals
  int<lower=0> K; // total number of covariates across all channels
  array[T] matrix[N, K] X; // covariates as an array of N x K matrices
  row_vector[K] X_m; // Means of all covariates at first time point
  int<lower=1> D; // number of B-splines
  matrix[D, T] Bs; // B-spline basis matrix
  int<lower=0> M; // number group-level effects (including intercepts)
  // number of fixed, varying and random coefficients, and related indices
  int<lower=0> K_varying_Hedelmallisyys;
  int<lower=0> K_random_Hedelmallisyys; // includes the potential random intercept
  int<lower=0> K_Hedelmallisyys; // K_fixed + K_varying
  array[K_varying_Hedelmallisyys] int J_varying_Hedelmallisyys;
  array[K_Hedelmallisyys] int J_Hedelmallisyys; // fixed and varying
  array[K_varying_Hedelmallisyys] int L_varying_Hedelmallisyys;
  // Parameters of vectorized priors
  matrix[K_varying_Hedelmallisyys, 2] delta_prior_pars_Hedelmallisyys;
  matrix[K_varying_Hedelmallisyys, 2] tau_prior_pars_Hedelmallisyys;
  matrix[K_random_Hedelmallisyys, 2] sigma_nu_prior_pars_Hedelmallisyys;
  matrix[N, T] y_Hedelmallisyys;
}
transformed data {
}
parameters {
  // Random group-level effects
  vector<lower=0>[M] sigma_nu; // standard deviations of random effects
  matrix[N, M] nu_raw;
  matrix[K_varying_Hedelmallisyys, D] omega_Hedelmallisyys; // Spline coefficients
  vector<lower=0.01>[K_varying_Hedelmallisyys] tau_Hedelmallisyys; // SDs for the random walks
  real a_Hedelmallisyys; // Mean of the first time point
  row_vector[D - 1] omega_raw_alpha_Hedelmallisyys; // Coefficients for alpha
  real<lower=0.01> tau_alpha_Hedelmallisyys; // SD for the random walk
  vector<lower=0>[N] sigma_Hedelmallisyys; // SD of the normal distribution
}
transformed parameters {
  vector[1] sigma_nu_Hedelmallisyys = sigma_nu[1:1];
  matrix[N, 1] nu_Hedelmallisyys = diag_post_multiply(nu_raw[, 1:1], sigma_nu_Hedelmallisyys);
  // Time-varying coefficients
  array[T] vector[K_varying_Hedelmallisyys] delta_Hedelmallisyys;
  // Time-varying intercept
  array[T] real alpha_Hedelmallisyys;
  // Spline coefficients
  real omega_alpha_1_Hedelmallisyys;
  row_vector[D] omega_alpha_Hedelmallisyys;
  for (t in 1:T) {
    delta_Hedelmallisyys[t] = omega_Hedelmallisyys * Bs[, t];
  }
  // Define the first alpha using mean a_Hedelmallisyys
  {
    vector[K_Hedelmallisyys] gamma__Hedelmallisyys;
    gamma__Hedelmallisyys[L_varying_Hedelmallisyys] = delta_Hedelmallisyys[1];
    omega_alpha_1_Hedelmallisyys = a_Hedelmallisyys - X_m[J_Hedelmallisyys] * gamma__Hedelmallisyys;
  }
  omega_alpha_Hedelmallisyys[1] = omega_alpha_1_Hedelmallisyys;
  omega_alpha_Hedelmallisyys[2:D] = omega_raw_alpha_Hedelmallisyys;
  for (t in 1:T) {
    alpha_Hedelmallisyys[t] = omega_alpha_Hedelmallisyys * Bs[, t];
  }
}
model {
  to_vector(nu_raw) ~ std_normal();
  sigma_nu_Hedelmallisyys ~ normal(sigma_nu_prior_pars_Hedelmallisyys[, 1], sigma_nu_prior_pars_Hedelmallisyys[, 2]);
  a_Hedelmallisyys ~ normal(2, 2);
  omega_raw_alpha_Hedelmallisyys[1] ~ normal(omega_alpha_1_Hedelmallisyys, tau_alpha_Hedelmallisyys);
  for (i in 2:(D - 1)) {
    omega_raw_alpha_Hedelmallisyys[i] ~ normal(omega_raw_alpha_Hedelmallisyys[i - 1], tau_alpha_Hedelmallisyys);
  }
  tau_alpha_Hedelmallisyys ~ normal(0, 2);
  omega_Hedelmallisyys[, 1] ~ normal(delta_prior_pars_Hedelmallisyys[, 1], delta_prior_pars_Hedelmallisyys[, 2]);
  for (i in 2:D) {
    omega_Hedelmallisyys[, i] ~ normal(omega_Hedelmallisyys[, i - 1], tau_Hedelmallisyys);
  }
  tau_Hedelmallisyys ~ normal(tau_prior_pars_Hedelmallisyys[, 1], tau_prior_pars_Hedelmallisyys[, 2]);
  sigma_Hedelmallisyys ~ exponential(1);
  {
    real ll = 0.0;
    vector[K_Hedelmallisyys] gamma__Hedelmallisyys;
    for (t in 1:T) {
      vector[N] intercept_Hedelmallisyys = alpha_Hedelmallisyys[t] + nu_Hedelmallisyys[, 1];
      gamma__Hedelmallisyys[L_varying_Hedelmallisyys] = delta_Hedelmallisyys[t];
      ll += normal_id_glm_lupdf(y_Hedelmallisyys[, t] | X[t][, J_Hedelmallisyys], intercept_Hedelmallisyys, gamma__Hedelmallisyys, sigma_Hedelmallisyys);
    }
    target += ll;
  }
}
generated quantities {
}
"
fit_malli_norm_var_its <- dynamite(
  dformula = malli_norm,
  data = ei_nollia_kunnat_skaalattu, 
  custom_stan_model = uusi_stan_koodi_norm_its,
  time = "Vuosi", group = "Kunta",
  chains = 4, cores = 4, seed = 189, refresh = 1, iter = 30000, warmup = 10000, thin=10
)

#Tallennetaan malli
#saveRDS(fit_malli_norm_var_its, "fit_malli_norm_var_its_lag.rds")
fit_malli_norm_var_its <- readRDS("fit_malli_norm_var_its_lag.rds")

#Konvergenssin tarkastelu
ess_r(fit_malli_norm_var_its)

#Sovitteen posteriorinäytteet
sovitteet_norm_var_its <- fitted(fit_malli_norm_var_its, thin = 4)

#Tallennetaan sovitteet
#saveRDS(sovitteet_norm_var_its, "sovitteet_norm_var_its.rds")
#sovitteet_norm_var_its <- readRDS("sovitteet_norm_var_its.rds")

#R2
norm_var_r2_its <-R2(sovitteet_norm_var_its)

#Posteriorikeskiarvo
mean(norm_var_r2_its)

#95% posterioriväli
quantile(norm_var_r2_its, probs = c(0.025, 0.975))

#--------------
#Gamma-mallit
#--------------
#
#Yhteinen phi
#
malli_gamma<-obs(Hedelmallisyys ~ -1 + varying(~ 1 + lag(Hedelmallisyys, k=1) + Maakunta + Korkea_aste_osuus + Perusaste_osuus + Tyottomat_osuus + Miehet_muut_tyon_ulkopuolella_osuus + Opiskelijat_osuus + Naiset_muut_tyon_ulkopuolella_osuus + Kunnan_ulkopuolella_tyossa + log_Vakiluku + Ei_uskontokuntaa + Keski_ika + Asuinalueella_syntyneet + Naimisissa_osuus + Yhden_henkilon_asuntokunnat_osuus + Neljan_henkilon_asuntokunnat_osuus + Vahintaan_seitseman_henkilon_asuntokunnat_osuus + Naiset_osuus_19_49) + random(~ 1),
                family = "gamma") + splines(df = 10, lb_tau = 0.01)

fit_malli_gamma <- dynamite(
  dformula = malli_gamma,
  data = ei_nollia_kunnat_skaalattu, time = "Vuosi", group = "Kunta",
  chains = 4, cores = 4, seed = 189, refresh = 1, iter = 30000, warmup = 10000, thin = 10
)

#Tallennetaan malli
#saveRDS(fit_malli_gamma, "fit_malli_gamma_lag.rds")
fit_malli_gamma <- readRDS("fit_malli_gamma_lag.rds")

#Konvergenssin tarkastelu
ess_r(fit_malli_gamma)

#Sovitteen posteriorinäytteet
sovitteet_gamma <- fitted(fit_malli_gamma, expand = TRUE, thin = 4)

#Tallennetaan sovitteet
#saveRDS(sovitteet_gamma, "sovitteet_gamma.rds")
#sovitteet_gamma <- readRDS("sovitteet_gamma.rds")

#R2
gamma_r2 <-R2(sovitteet_gamma)

#Posteriorikeskiarvo
mean(gamma_r2)

#95% posterioriväli
quantile(gamma_r2, probs = c(0.025, 0.975))


#
#Kuntakohtainen phi
#
uusi_stan_koodi_gamma <- "
functions {
}
data {
  int<lower=1> T; // number of time points
  int<lower=1> N; // number of individuals
  int<lower=0> K; // total number of covariates across all channels
  array[T] matrix[N, K] X; // covariates as an array of N x K matrices
  row_vector[K] X_m; // Means of all covariates at first time point
  int<lower=1> D; // number of B-splines
  matrix[D, T] Bs; // B-spline basis matrix
  int<lower=0> M; // number group-level effects (including intercepts)
  // number of fixed, varying and random coefficients, and related indices
  int<lower=0> K_varying_Hedelmallisyys;
  int<lower=0> K_random_Hedelmallisyys; // includes the potential random intercept
  int<lower=0> K_Hedelmallisyys; // K_fixed + K_varying
  array[K_varying_Hedelmallisyys] int J_varying_Hedelmallisyys;
  array[K_Hedelmallisyys] int J_Hedelmallisyys; // fixed and varying
  array[K_varying_Hedelmallisyys] int L_varying_Hedelmallisyys;
  // Parameters of vectorized priors
  matrix[K_varying_Hedelmallisyys, 2] delta_prior_pars_Hedelmallisyys;
  matrix[K_varying_Hedelmallisyys, 2] tau_prior_pars_Hedelmallisyys;
  matrix[K_random_Hedelmallisyys, 2] sigma_nu_prior_pars_Hedelmallisyys;
  matrix<lower=0>[N, T] y_Hedelmallisyys;
}
transformed data {
}
parameters {
  // Random group-level effects
  vector<lower=0>[M] sigma_nu; // standard deviations of random effects
  matrix[N, M] nu_raw;
  matrix[K_varying_Hedelmallisyys, D] omega_Hedelmallisyys; // Spline coefficients
  vector<lower=0.01>[K_varying_Hedelmallisyys] tau_Hedelmallisyys; // SDs for the random walks
  real a_Hedelmallisyys; // Mean of the first time point
  row_vector[D - 1] omega_raw_alpha_Hedelmallisyys; // Coefficients for alpha
  real<lower=0.01> tau_alpha_Hedelmallisyys; // SD for the random walk
  vector<lower=0>[N] phi_Hedelmallisyys; // Shape parameter of the Gamma distribution
}
transformed parameters {
  vector[1] sigma_nu_Hedelmallisyys = sigma_nu[1:1];
  matrix[N, 1] nu_Hedelmallisyys = diag_post_multiply(nu_raw[, 1:1], sigma_nu_Hedelmallisyys);
  // Time-varying coefficients
  array[T] vector[K_varying_Hedelmallisyys] delta_Hedelmallisyys;
  // Time-varying intercept
  array[T] real alpha_Hedelmallisyys;
  // Spline coefficients
  real omega_alpha_1_Hedelmallisyys;
  row_vector[D] omega_alpha_Hedelmallisyys;
  for (t in 1:T) {
    delta_Hedelmallisyys[t] = omega_Hedelmallisyys * Bs[, t];
  }
  // Define the first alpha using mean a_Hedelmallisyys
  {
    vector[K_Hedelmallisyys] gamma__Hedelmallisyys;
    gamma__Hedelmallisyys[L_varying_Hedelmallisyys] = delta_Hedelmallisyys[1];
    omega_alpha_1_Hedelmallisyys = a_Hedelmallisyys - X_m[J_Hedelmallisyys] * gamma__Hedelmallisyys;
  }
  omega_alpha_Hedelmallisyys[1] = omega_alpha_1_Hedelmallisyys;
  omega_alpha_Hedelmallisyys[2:D] = omega_raw_alpha_Hedelmallisyys;
  for (t in 1:T) {
    alpha_Hedelmallisyys[t] = omega_alpha_Hedelmallisyys * Bs[, t];
  }
}
model {
  to_vector(nu_raw) ~ std_normal();
  sigma_nu_Hedelmallisyys ~ normal(sigma_nu_prior_pars_Hedelmallisyys[, 1], sigma_nu_prior_pars_Hedelmallisyys[, 2]);
  a_Hedelmallisyys ~ normal(0.69, 2);
  omega_raw_alpha_Hedelmallisyys[1] ~ normal(omega_alpha_1_Hedelmallisyys, tau_alpha_Hedelmallisyys);
  for (i in 2:(D - 1)) {
    omega_raw_alpha_Hedelmallisyys[i] ~ normal(omega_raw_alpha_Hedelmallisyys[i - 1], tau_alpha_Hedelmallisyys);
  }
  tau_alpha_Hedelmallisyys ~ normal(0, 2);
  omega_Hedelmallisyys[, 1] ~ normal(delta_prior_pars_Hedelmallisyys[, 1], delta_prior_pars_Hedelmallisyys[, 2]);
  for (i in 2:D) {
    omega_Hedelmallisyys[, i] ~ normal(omega_Hedelmallisyys[, i - 1], tau_Hedelmallisyys);
  }
  tau_Hedelmallisyys ~ normal(tau_prior_pars_Hedelmallisyys[, 1], tau_prior_pars_Hedelmallisyys[, 2]);
  phi_Hedelmallisyys ~ exponential(1);
  {
    real ll = 0.0;
    for (t in 1:T) {
      vector[N] intercept_Hedelmallisyys = alpha_Hedelmallisyys[t] + nu_Hedelmallisyys[, 1] + X[t][, J_varying_Hedelmallisyys] * delta_Hedelmallisyys[t];
      ll += gamma_lupdf(y_Hedelmallisyys[, t] | phi_Hedelmallisyys, phi_Hedelmallisyys .* inv(exp(intercept_Hedelmallisyys)));
    }
    target += ll;
  }
}
generated quantities {
}
"
fit_malli_gamma_var <- dynamite(
  dformula = malli_gamma,
  custom_stan_model = uusi_stan_koodi_gamma,
  data = ei_nollia_kunnat_skaalattu, time = "Vuosi", group = "Kunta",
  chains = 4, cores = 4, seed = 189, refresh = 1, iter = 30000, warmup = 10000, thin = 10
)

#Tallennetaan malli
#saveRDS(fit_malli_gamma_var, "fit_malli_gamma_var_lag.rds")
fit_malli_gamma_var <- readRDS("fit_malli_gamma_var_lag.rds")

#Konvergenssin tarkastelu
ess_r(fit_malli_gamma_var)

#Sovitteen posteriorinäytteet
sovitteet_gamma_var <- fitted(fit_malli_gamma_var, expand = TRUE, thin = 4)

#Tallennetaan sovitteet
#saveRDS(sovitteet_gamma_var, "sovitteet_gamma_var.rds")
#sovitteet_gamma_var <- readRDS("sovitteet_gamma_var.rds")

#R2
gamma_var_r2 <-R2(sovitteet_gamma_var)

#Posteriorikeskiarvo
mean(gamma_var_r2)

#95% posterioriväli
quantile(gamma_var_r2, probs = c(0.025, 0.975))

#--------------
#Yksi-pois-ristiinvalidointi
#--------------

#Toteutetaan yksi-pois-ristiinvalidointi
set.seed(189)
loo_norm<-loo(fit_malli_norm, thin = 4)
loo_norm_var<-loo(fit_malli_norm_var, thin = 4)
loo_norm_var_its<-loo(fit_malli_norm_var_its, thin = 4)
loo_gamma<-loo(fit_malli_gamma, thin = 4)
loo_gamma_var<-loo(fit_malli_gamma_var, thin = 4)

#Tulosten tallennus ja lataus
#saveRDS(loo_norm, "loo_norm.rds")
#loo_norm <- readRDS("loo_norm.rds")

#saveRDS(loo_norm_var, "loo_norm_var.rds")
#loo_norm_var <- readRDS("loo_norm_var.rds")

#saveRDS(loo_norm_var_its, "loo_norm_var_its.rds")
#loo_norm_var_its <- readRDS("loo_norm_var_its.rds")

#saveRDS(loo_gamma, "loo_gamma.rds")
#loo_gamma <- readRDS("loo_gamma.rds")

#saveRDS(loo_gamma_var, "loo_gamma_var.rds")
#loo_gamma_var <- readRDS("loo_gamma_var.rds")

#Vertaillaan tuloksia
vert <- loo::loo_compare(x = list(loo_norm, loo_norm_var, loo_norm_var_its, loo_gamma, loo_gamma_var))
print(vert, simplify = FALSE)