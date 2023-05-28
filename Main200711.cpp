/**���ÿ�**/
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

/**�����ܻ���**/
double Patm = 1.01325e5;
double hfg = 2500;
double InterP(int index, double A[], double B[], double a);                                                      /**��׼����ѹ**/

int main(){

    /**����\����ļ�-���������������**/
    ifstream WRrelation("0-1-Isotherm-MOF-303.txt");                  /**���ϵ�����������**/
    ifstream PTrelation("0-2-SaturatedPressure.txt");                /**�¶�-����ѹ����**/
    ofstream CErelation;
    CErelation.open("1-output.txt");

    ofstream Output;
    Output.open("1-FinalOutput.txt");
    ofstream Output0, Output1, Output2, Output3;
    Output0.open("2-FinalOutput0.txt");
    Output1.open("2-FinalOutput1.txt");
    Output2.open("2-FinalOutput2.txt");
    Output3.open("2-FinalOutput3.txt");

    int WRnum, PTnum, p;
    WRrelation>>WRnum;
    PTrelation>>PTnum;
    double MWDa[WRnum], MRHa[WRnum];
    double MWDd[WRnum], MRHd[WRnum];
    double MT[PTnum], MPvs[PTnum];
    for(p = 0; p <= WRnum; p++){WRrelation>>MRHa[p]>>MWDa[p];if(p <= WRnum)WRrelation>>MRHd[p]>>MWDd[p];}
    for(p = 0; p <= PTnum; p++){PTrelation>>MT[p]>>MPvs[p];}

    /**���������**/
    /**��װ�������й�**/
    double WDcycle = 0.2;                                    /**������ѭ��������**/
    double MassRatio = 1, Qst = 2889, cps = 1;               /**����������**/
    double FlowRatio;

    /**��װ�ò�������й�**/
    //double DLTG_ADout = 5;                                 /**���������ڿ���������Դ���¶�**/
    double DLTS_ADend = 5;                                   /**����������ĩ�ڸ�����Դ���¶�**/
    double DLYS_ADend = 1;                                   /**����������ĩ�ڵ�����ڿ����ľ���ʪ�ȣ������Ʋ**/
    double DLTS_DE = 5;                                      /**��������������������������Դ���¶ȣ������Ʋ**/
    double DLYS_DE = 2;                                      /**���������������������׵��ڸ��ڿ������������ľ���ʪ��**/
    double DLTG_CDout = 5;                                   /**���������ڿ���������Դ���¶ȣ������Ʋ**/
    //double DLTG_DEout = 5;                                 /**�������������̳��ڿ����¶ȵ�����Դ���¶�**/
    //double TytS_AD = 0.95;                                   /**�����������¶�����ڿ������µ�Ӱ��**/
    double TytG_AD = 0.75;                                   /**�������������ڿ����¶����������¶�Ӱ��**/
    double TytG_DE = 0.80;                                   /**�������������ڿ����¶����������¶�Ӱ��**/
    //double Yyt_AD = 0.1;

    /**�����ܶ�**/
    int TGin_MIN = 0, TGin_MAX = 50, TGin_DLT = 2;               /**�������״̬���¶Ƚڵ�**/
    double YGrate_MIN = 0.5, YGrate_MAX = 25, YGrate_DLT = 0.01;
    double TL_AD_MIN = -5, TL_AD_DLT = 1;
    double TL_CD_MIN = -5, TL_CD_DLT = 1;
    int WDnum = 100;

    /**���������м���**/
    double YGrate, RHG_ADin;                                  /**ȡˮ����**/
    int TG_ADin;
    double TL_AD;
    double TL_CD;

    double WD_ADend, WD_DEend, WD_DEtemp;
    double RHS_ADend, RHS_DEend, RHS_DEavrg;
    double TG_ADout, TG_CDout, TS_ADend, TS_DE, TS_DEcheck, TL_DE, TG_DEout;
    double PG_CDout, PS_ADend, PS_DEavrg, PS_DEend, PG_ADin;
    double YG_ADin, YG_CDout, YS_ADend, YS_ADavrg, YS_DEavrg, YG_DEout;
    double DLI1, DLI2, DLI3;
    double Exergy;

    double TL_CD_OP, TL_AD_OP, TL_DE_OP, Exergy_OP, YGrate_OP, FlowRatio_OP;
    double TL_CD_OP1, TL_AD_OP1, TL_DE_OP1, Exergy_OP1, YGrate_OP1, FlowRatio_OP1;
    double TL_CD_OP2, TL_DE_OP2, Exergy_OP2;

    for (RHG_ADin = 0; RHG_ADin <= 100; RHG_ADin = RHG_ADin + 0.1){

        cout<<RHG_ADin<<endl;
        CErelation<<RHG_ADin;
        Output<<RHG_ADin;
        Output0<<RHG_ADin;
        Output1<<RHG_ADin;
        Output2<<RHG_ADin;
        Output3<<RHG_ADin;

        for (TG_ADin = TGin_MIN; TG_ADin <= TGin_MAX + 1E-3; TG_ADin = TG_ADin + TGin_DLT){

            /**��������**/
            PG_ADin = exp(-5800.2206 / (TG_ADin + 273.15) + 1.3914993 - 0.048640239 * (TG_ADin + 273.15) + 0.000041764768 * pow(TG_ADin + 273.15,2) - 0.000000014452093 * pow(TG_ADin + 273.15,3) + 6.5459673 * log(TG_ADin + 273.15));
                                                                                                   /**������AD��ڿ���G��������ѹ**/
            YG_ADin =  622 * RHG_ADin / 100 / (Patm / PG_ADin - RHG_ADin / 100);                   /**������AD��ڿ���Gʪ��**/

            /**�������������Դ�¶�**/
            Exergy_OP = -10;
            for (TL_AD = TL_AD_MIN; TL_AD <= TG_ADin + 1E-3; TL_AD = TL_AD + TL_AD_DLT){
            //for (TL_AD = TG_ADin; TL_AD <= TG_ADin + 1E-3; TL_AD = TL_AD + TL_AD_DLT){
                Exergy_OP1 = -10;
                for (YGrate = YGrate_MIN; YGrate <= YGrate_MAX + 1E-4; YGrate = YGrate + YGrate_DLT){
                //for (YGrate = 8; YGrate <= 8 + 1E-4; YGrate = YGrate + YGrate_DLT){

                    DLYS_DE = 0.4 * YGrate;
                    DLYS_ADend = 0.2 * YGrate;
                    DLTS_DE = 6 + 0.2 * YGrate;
                    DLTS_ADend = 4 + 0.2 * YGrate;
                    DLTG_CDout = 4 + 0.2 * YGrate;

                    /**�������̼���**/
                    TS_ADend = TL_AD + DLTS_ADend;
                    //TS_ADend = (TytS_AD * TL_AD + (1 - TytS_AD) * TG_ADin >= TS_ADend) ? (TytS_AD * TL_AD + (1 - TytS_AD) * TG_ADin) : TS_ADend;
                    //������ڿ������¶����������ȵ�Ӱ��//
                                                                                                        /**������S����ĩ���¶�**/
                    YS_ADend = YG_ADin - DLYS_ADend;                                                    /**������S����ĩ�ڿ׵���ʪ��**/
                    PS_ADend = exp(-5800.2206 / (TS_ADend + 273.15) + 1.3914993 - 0.048640239 * (TS_ADend + 273.15) + 0.000041764768 * pow(TS_ADend + 273.15,2) - 0.000000014452093 * pow(TS_ADend + 273.15,3) + 6.5459673 * log(TS_ADend + 273.15));
                                                                                                        /**������S����ĩ�ڿ׵����ͺ�ʪ��**/
                    RHS_ADend = Patm * 100 / (622 * PS_ADend / YS_ADend + PS_ADend);                    /**������S����ĩ�ڿ׵����ʪ��**/
                    if (RHS_ADend <= 0) {continue;}
                    WD_ADend = InterP(WRnum, MRHa, MWDa, RHS_ADend);                                    /**������S����ĩ��������**/

                    /**��������**/
                    WD_DEend = WD_ADend - WDcycle;                                                      /**������S����ĩ��������**/
                    if (WD_DEend <= 0) {continue;}
                    RHS_DEend = InterP(WRnum, MWDd, MRHd, WD_DEend);                                    /**������S����ĩ�ڿ׵����ʪ��**/
                    RHS_DEavrg = 0, WD_DEtemp = WD_DEend;
                    for (p = 0; p < WDnum; p++){RHS_DEavrg = InterP(WRnum, MWDd, MRHd, WD_DEtemp) / WDnum + RHS_DEavrg, WD_DEtemp = WD_DEtemp + WDcycle / WDnum;}
                                                                                                    /**������S�������̿׵�ƽ�����ʪ��**/
                    //cout<<RHG_ADin<<'\t'<<WD_ADend<<'\t'<<RHS_DEend<<'\t'<<RHS_DEavrg<<endl;
                    YS_ADavrg =  622 * RHS_DEavrg / (Patm / PS_ADend - RHS_DEavrg / 100) / 100;
                    FlowRatio = max((YGrate + DLYS_DE) / (YG_ADin - YS_ADavrg),2.0);
                    TytG_AD = 0.8 / pow(FlowRatio,0.2);
                    TG_ADout = (TG_ADin <= TS_ADend) ? TG_ADin : TytG_AD * TS_ADend + (1 - TytG_AD) * TG_ADin;
                    //cout<<TG_ADin<<'\t'<<YG_ADin<<'\t'<<YS_ADend<<'\t'<<RHS_ADend<<'\t'<<RHS_DEavrg<<'\t'<<YS_ADavrg<<'\t'<<FlowRatio<<'\t'<<TytG_AD<<endl;

                    /**��������**/
                    Exergy_OP2 = -10;
                    for (TL_CD = TL_CD_MIN; TL_CD <= TG_ADin + 1E-3; TL_CD = TL_CD + TL_CD_DLT){
                    //for (TL_CD = TG_ADin; TL_CD <= TG_ADin + 1E-3; TL_CD = TL_CD + TL_CD_DLT){

                        TG_CDout = TL_CD + DLTG_CDout;                                                    /**���������ڿ���G�¶�**/
                        //cout<<TL_CD<<endl;
                        PG_CDout = exp(-5800.2206 / (TG_CDout + 273.15) + 1.3914993 - 0.048640239 * (TG_CDout + 273.15) + 0.000041764768 * pow(TG_CDout + 273.15,2) - 0.000000014452093 * pow(TG_CDout + 273.15,3) + 6.5459673 * log(TG_CDout + 273.15));
                                                                                                          /**���������ڿ���G��������ѹ**/
                        YG_CDout =  622 / (Patm / PG_CDout - 1);                                            /**���������ڿ���G��ʪ��**/

                        /**��������**/
                        YG_DEout = YG_CDout + YGrate;                                                       /**�������̳��ڿ���Gƽ����ʪ��**/
                        YS_DEavrg = YG_DEout + DLYS_DE;                                                     /**��������������S�׵���ƽ����ʪ��**/
                        PS_DEavrg = Patm * 100 / (622 * RHS_DEavrg / YS_DEavrg + RHS_DEavrg);               /**��������������S�׵��ڱ�������ѹ**/
                        TS_DE = InterP(PTnum, MPvs, MT, PS_DEavrg / 1000);                                  /**��������������S�¶�**/

                        PS_DEend = Patm * 100 / (622 * RHS_DEend / YG_CDout + RHS_DEend);                   /**����ĩ��������S�׵��ڱ�������ѹ**/
                        TS_DEcheck = InterP(PTnum, MPvs, MT, PS_DEend / 1000);                              /**����ĩ��������S�¶�**/
                        //У�������������ܷ񽫻�����������RHS_DEend//

                        TL_DE = max(TS_DE, TS_DEcheck) + DLTS_DE;                                           /**����������ԴS�¶�**/
                        if (TL_DE >= 100) {continue;}
                        TG_DEout = (TL_DE - DLTS_DE >= TG_CDout) ? TytG_DE * (TL_DE - DLTS_DE) + (1 - TytG_DE) * TG_CDout : TG_CDout;         /**���������ڿ���G�¶�**/

                        //cout<<"attention"<<'\t'<<TG_ADin<<'\t'<<TG_CDout<<'\t'<<PG_CDout<<'\t'<<Ydesic<<'\t'<<PS_DEavrg<<'\t'<<Tdesic<<endl;

                        DLI1 = TG_DEout * 1.01 + YG_DEout / 1000 * (Qst + 1.84 * TG_DEout) - TG_CDout * 1.01 - YG_CDout / 1000 * (Qst + 1.84 * TG_CDout) + YGrate / WDcycle / 1000 * (TS_DE - TS_ADend) * cps + YGrate / WDcycle / 1000 * (TL_DE - TL_AD) * cps * MassRatio;


                        DLI2 = FlowRatio * (TG_ADin * 1.01 + YG_ADin / 1000 * (Qst + 1.84 * TG_ADin) - TG_ADout * 1.01 - (YG_ADin - YGrate / FlowRatio) / 1000 * (Qst + 1.84 * TG_ADout)) + YGrate / WDcycle / 1000 * (TS_DE - TS_ADend) * cps + YGrate / WDcycle / 1000 * (TL_DE - TL_AD) * cps * MassRatio;
                        DLI3 = TG_DEout * 1.01 + YG_DEout / 1000 * (hfg + 1.84 * TG_DEout) - TG_CDout * 1.01 - YG_CDout / 1000 * (hfg + 1.84 * TG_CDout);
                        Exergy = (DLI1 * (TL_DE - TG_ADin) / (TL_DE + 273.15) + DLI2 * (TG_ADin - TL_AD) / (TL_AD + 273.15) + DLI3 * (TG_ADin - TL_CD) / (TL_CD + 273.15)) * 1000 / YGrate;

                        //if(RHG_ADin == 40 && TG_ADin == 10)cout<<TL_AD<<'\t'<<PG_ADin<<'\t'<<YG_ADin<<'\t'<<TS_ADend<<'\t'<<YS_ADend<<'\t'<<PS_ADend<<'\t'<<RHS_ADend<<'\t'<<WD_ADend<<endl;
                        //if(RHG_ADin == 40 && TG_ADin == 10)cout<<WD_ADend<<'\t'<<TG_ADout<<'\t'<<InterP(WRnum, MWDd, MRHd, WD_ADend)<<'\t'<<WD_DEend<<'\t'<<RHS_DEend<<'\t'<<RHS_DEavrg<<'\t'<<PG_CDout<<'\t'<<YG_CDout<<endl;
                        //if(RHG_ADin == 46 && TG_ADin == 6)cout<<YG_CDout<<'\t'<<YG_DEout<<'\t'<<YS_DEavrg<<'\t'<<RHS_DEavrg<<'\t'<<PS_DEavrg<<'\t'<<TS_DE<<'\t'<<PS_DEend<<'\t'<<TS_DEcheck<<'\t'<<TL_DE<<endl;
                        //cout<<TG_ADin<<'\t'<<TL_AD<<'\t'<<TL_CD<<'\t'<<TL_DE<<'\t'<<TG_DEout<<'\t'<<DLI1<<'\t'<<DLI2<<'\t'<<DLI3<<'\t'<<DLI1 * (TL_DE - TG_ADin) / (TL_DE + 273.15)<<'\t'<<DLI2 * (TG_ADin - TL_AD) / (TL_AD + 273.15)<<'\t'<<DLI2 * (TG_ADin - TL_CD) / (TL_CD + 273.15)<<'\t'<<Exergy<<endl;

                        if (Exergy_OP2 < 0 || Exergy <= Exergy_OP2){Exergy_OP2 = Exergy, TL_CD_OP2 = TL_CD, TL_DE_OP2 = TL_DE;}
                    }
                    if (Exergy_OP2 > 0 && (Exergy_OP1 < 0 || Exergy_OP2 <= Exergy_OP1)){Exergy_OP1 = Exergy_OP2, TL_CD_OP1 = TL_CD_OP2, TL_DE_OP1 = TL_DE_OP2, YGrate_OP1 = YGrate; FlowRatio_OP1 = FlowRatio;}
                }
                if (Exergy_OP1 > 0 && (Exergy_OP < 0 || Exergy_OP1 <= Exergy_OP)){Exergy_OP = Exergy_OP1, TL_CD_OP = TL_CD_OP1, TL_DE_OP = TL_DE_OP1, YGrate_OP = YGrate_OP1, FlowRatio_OP = FlowRatio_OP1, Exergy_OP = Exergy_OP1, TL_AD_OP = TL_AD;}
            }

            //if (Exergy_OP < 0 || TL_DE_OP >= 95){
            if (Exergy_OP < 0){
                Output<<'\t'<<0.0/0.0<<'\t'<<0.0/0.0;
                Output0<<'\t'<<0.0/0.0<<'\t'<<0.0/0.0;
                Output1<<'\t'<<0.0/0.0<<'\t'<<0.0/0.0;
                Output2<<'\t'<<0.0/0.0<<'\t'<<0.0/0.0;
                Output3<<'\t'<<0.0/0.0<<'\t'<<0.0/0.0;
                CErelation<<'\t'<<TG_ADin<<'\t'<<0.0/0.0<<'\t'<<0.0/0.0<<'\t'<<0.0/0.0<<'\t'<<0.0/0.0<<'\t'<<0.0/0.0;
            }
            else {
                Output<<'\t'<<YG_ADin / 1000<<'\t'<<Exergy_OP;
                Output0<<'\t'<<YG_ADin / 1000<<'\t'<<YGrate_OP;
                Output1<<'\t'<<YG_ADin / 1000<<'\t'<<TL_AD_OP;
                Output2<<'\t'<<YG_ADin / 1000<<'\t'<<TL_CD_OP;
                Output3<<'\t'<<YG_ADin / 1000<<'\t'<<TL_DE_OP;
                CErelation<<'\t'<<TG_ADin<<'\t'<<TL_AD_OP<<'\t'<<TL_CD_OP<<'\t'<<TL_DE_OP<<'\t'<<YGrate_OP<<'\t'<<Exergy_OP<<'\t'<<FlowRatio;
            }
        }
        CErelation<<endl;
        Output<<endl;
        Output0<<endl;
        Output1<<endl;
        Output2<<endl;
        Output3<<endl;
    }

    WRrelation.close();
    PTrelation.close();
    CErelation.close();
    Output.close();
    Output0.close();
    Output1.close();
    Output2.close();
    Output3.close();
    return 0;
}

double InterP(int index, double A[], double B[], double a){

    double b;
    int i;
    if (a < A[0]) {b = 0.0 / 0.0;}
    else for(i = 0; i < index; i++){if(A[i] <= a && a <= A[i + 1])break;}
    if(i == index){b = B[index];}
    else b = (a - A[i]) / (A[i + 1] - A[i]) * (B[i + 1] - B[i]) + B[i];
    return b;
}
