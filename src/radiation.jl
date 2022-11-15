# alb = similar(albs)   ### hack
# asrf_out = similar(albs)   ### hack
# Sdirt = similar(albs)   ### hack
# Sdift = similar(albs)   ### hack
# SWveg = similar(albs)   ### hack
# SWsrf = similar(albs)   ### hack
# SWsci = similar(albs)   ### hack
# LWt = similar(albs)   ### hack
# SWtopo_out = similar(albs)   ### hack

# Snow albedo
for j = 1:Ny
    for i = 1:Nx

        if (tilefrac[i,j] >= tthresh) # exclude points outside tile of interest

            # New Snow albedo 
            if (OSHDTN == 0)
                global afs ### hack
                afs = asmx
            else # OSHDTN == 1
                # 11/2021 tuning: high elevation afs changed from 0.86 to 0.92
                if (dem[i,j] >= 2200)
                    afs  = 0.92
                elseif (dem[i,j] <= 1400)
                    afs = 0.80
                else
                    afs = 0.92 + (2200 - dem[i,j]) / (2200 - 1400) * (0.80 - 0.92)
                end
            end

            if (ALBEDO == 0)
                # Diagnostic
                albs[i,j] = asmn + (afs - asmn)*(Tsrf[i,j] - Tm) / Talb
                if (albs[i,j] < min(afs, asmn))
                    albs[i,j] = min(afs, asmn)
                end
                if (albs[i,j] > max(afs, asmn))
                    albs[i,j] = max(afs, asmn)
                end

            elseif (ALBEDO == 1)
                # Prognostic
                tau = tcld
                if (Tsrf[i,j] >= Tm)
                    tau = tmlt
                end
                # Forest adjustments -> not yet properly tested for OSHD but option currently unused
                if (month > 4 && month < 10)
                  tau = 70*3600
                end

                if fveg[i,j] > 0 && Sdir[i,j] > eps(Float64)
                  tau = tau /((1-trcn[i,j]*fsky[i,j])*(1+adfl*Tv[i,j]) + adfs*Tv[i,j])
                elseif fveg[i,j] > 0 && Sdif[i,j] > eps(Float64)
                  tau = tau/((1-trcn[i,j]*fsky[i,j]) + adfs*trcn[i,j]*fsky[i,j])
                elseif (fveg[i,j] > 0 && (Sdir[i,j] + Sdif[i,j] <= eps(Float64)))
                  tau = tau/(2-trcn[i,j]*fsky[i,j])
                end   

                rt = 1/tau + Sf[i,j]/Sfmin
                alim = (asmn/tau + Sf[i,j]*afs/Sfmin)/rt
                albs[i,j] = alim + (albs[i,j] - alim)*exp(-rt*dt)
                if (albs[i,j] < min(afs, asmn))
                    albs[i,j] = min(afs, asmn)
                end
                if (albs[i,j] > max(afs, asmn))
                    albs[i,j] = max(afs, asmn)
                end

            else # ALBEDO == 2
                # Prognostic, tuned, copied from JIM
                SWEtmp = sum(Sice[:,i,j] + Sliq[:,i,j])
                ## DECAY RATES
                # Melting and cold snow decay times
                if (OSHDTN==0)
                    adm=100
                    adc = 1000
                else
                    if (month > 6 && month < 10)
                        adm = 50
                    else
                        adm = 130
                    end
                    adc = 3000
                end
                if (Tsrf[i,j] >= Tm && Tsnow[1,i,j] >= Tm)  # was only based on Tss
                  albs[i,j] = (albs[i,j] - asmn)*exp(-(dt/3600)/adm) + asmn
                else
                  albs[i,j] = albs[i,j] - (dt/3600)/adc
                end
                if (SWEtmp < 75.0) # more stuff showing on and up through snow
                  afs = afs * 0.80
                end
                # Reset to fresh snow albedo (wasn't originally available; only else term)
                if ((Sf[i,j] * dt) > 0.0 && Sf24h[i,j] > Sfmin)
                  albs[i,j] = afs
                else
                  albs[i,j] = albs[i,j] + (afs - albs[i,j])*Sf[i,j]*dt/Sfmin
                end
                ## End Adjustments
                if (albs[i,j] > afs)
                    albs[i,j] = afs
                end
                if (albs[i,j] < asmn)
                    albs[i,j] = asmn
                end
              
            end
        end
    end
end

# Surface and canopy net shortwave radiation
for j = 1:Ny
    for i = 1:Nx
        if (tilefrac[i,j] >= tthresh) # exclude points outside tile of interest

            # Surface albedo
            asrf = albs[i,j]*(1-fveg[i,j]*fsar)  
            if (fsnow[i,j] <= eps(Float64))
                asrf = alb0[i,j]
            end

            # Partial snowcover on canopy
            fcans = 0
            if (scap[i,j] > eps(Float64))
                fcans = Sveg[i,j] / scap[i,j]
            end
            aveg = (1 - fcans)*avg0 + fcans*avgs
            acan = fveg[i,j]*aveg
            # Canopy surface albedo for computing terrain radiation over canopy
            alb[i,j] = fveg[i,j]*aveg + (1-fveg[i,j])*asrf
        
            # Surface albedo is stored in asurf_out to write in results
            asrf_out[i,j] = alb[i,j]
            
            #   hack not used
            #   if (RADSBG == 1)
            #     # Call Subgrid parameterization for SW radiation to compute SWtopo,netto SWtn
            #     call SWRADTOPO(alb[i,j],Sdir[i,j],Sdif[i,j],SWsrf[i,j],Sdirt[i,j],Sdift[i,j],SWtopo_out,Sun_elev,year,month,day,hour,i,j)
            #   end

            if (RADSBG == 0)
                global SWtopo_out  ### hack
                SWtopo_out =  alb[i,j]*(Sdir[i,j]+Sdif[i,j])
                Sdirt[i,j] = Sdir[i,j]
                Sdift[i,j] = Sdif[i,j] 
            end

            # Solar radiation trasmission 
            if (CANMOD == 0)
                SWveg[i,j] = 0
                SWsrf[i,j] = (1 - alb[i,j])*(Sdir[i,j]+Sdif[i,j])
                SWsci[i,j] = Sdift[i,j]+Sdir[i,j]  
            end

            if (CANMOD == 1)
                Sdif_aux = fsky[i,j]/fsky_terr[i,j]*Sdift[i,j]
                tdif = trcn[i,j]
                tdir = Tv[i,j] 

                # Effective albedo and net radiation
                alb[i,j] = acan + (1 - acan)*asrf*tdif^2
                if (Sdif_aux + Sdirt[i,j] > eps(Float64))
                alb[i,j] = (acan*(Sdif_aux+tdir*Sdirt[i,j]) + asrf*tdif*(tdif*Sdif_aux+tdir*Sdirt[i,j])) / (Sdif_aux + Sdirt[i,j])
                end
                SWsrf[i,j] = (1 - asrf)*(tdif*Sdif_aux + tdir*Sdirt[i,j])
                SWveg[i,j] = ((1-tdif)*(1-aveg)+tdif*asrf*(1-tdif))*Sdif_aux + (tdir*fveg[i,j]*(1-aveg)+tdir*asrf*(1-tdif))*Sdirt[i,j]   # local SWR absorption by vegetation correlates with local tdir  
                SWsci[i,j] = tdif*Sdif_aux + tdir*Sdirt[i,j]
            end 
        
            # Thermal emissions from surroundings 
            # Terrain LWR if not calculated later; 
            LWt[i,j] = fsky_terr[i,j]*LW[i,j] + (1 - fsky_terr[i,j])*sb*Ta[i,j]^4   

            # LW overwritten by LWt only if EBALFOR is used, where terrain impacts are accounted for already 
            if (CANMOD == 0 || fveg[i,j] == 0)
                LW[i,j] = LWt[i,j]
            end
                    
        end 
    
    end

end