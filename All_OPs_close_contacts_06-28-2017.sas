/* Multiple OP calculations*/
%let drive = E;
%let path = &drive.:\OneDrive\SAS_work\MD_SAS_analysis\All_plots\MultiPlots;
libname mltplt "&path.";

%let DTP_path = &drive.:\OneDrive\SAS_work\MD_SAS_analysis\DTP_analysis\working;
%let DCV_path = &drive.:\OneDrive\SAS_work\MD_SAS_analysis\DCV_analysis\working;
%let EFP_path = &drive.:\OneDrive\SAS_work\MD_SAS_analysis\EFP_analysis\working;
%let EPP_path = &drive.:\OneDrive\SAS_work\MD_SAS_analysis\EPP_analysis\working;
%let MAP_path = &drive.:\OneDrive\SAS_work\MD_SAS_analysis\MAP_analysis;
%let MVP_path = &drive.:\OneDrive\SAS_work\MD_SAS_analysis\MVP_analysis;
%let ODM_path = &drive.:\OneDrive\SAS_work\MD_SAS_analysis\ODM_analysis;
%let OMT_path = &drive.:\OneDrive\SAS_work\MD_SAS_analysis\OMT_analysis_test;
%let TMT_path = &drive.:\OneDrive\SAS_work\MD_SAS_analysis\TMT_analysis\working;
%let VMT_path = &drive.:\OneDrive\SAS_work\MD_SAS_analysis\VMT_analysis\working;

libname DTPlib "&DTP_path";
libname DCVlib "&DCV_path";
libname EFPlib "&EFP_path";
libname EPPlib "&EPP_path";
libname MAPlib "&MAP_path";
libname MVPlib "&MVP_path";
libname ODMlib "&ODM_path";
libname OMTlib "&OMT_path";
libname TMTlib "&TMT_path";
libname VMTlib "&VMT_path";


%let lignms = DTP DCV EFP EPP MAP MVP ODM OMT TMT VMT;

options minoperator orientation = landscape;
ods graphics on;
OPTIONS NODATE cpucount=3 threads;
/*===============================================================================*/

* This macro finds the length of a given delimiter separated list;
    %macro get_length(valst, sepr);
        %let valst = %qtrim(&valst);
        %if %length(valst) ne 0 %then %do;
            %let i = 1;
            %let varnm = %scan(&valst, &i, %str(&sepr));
            %do %until (&varnm=);
                %let varnm = %scan(&valst, &i, %str(&sepr));
                
                %let i = %eval(&i+1);
                %*put &i;               
            %end;
        %let i = %eval(&i - 2); %* to adjust;
        %end;
        %else
            %let i = 0;
        
        &i;
    %mend get_length;
/*===============================================================================*/
/*===============================================================================*/

%macro close_contact_vertic_concat(ligID, dst_id_lst, libry);
    %let no_of_files = %get_length(&dst_id_lst, %str( ));
    %let no_of_files = %qsysfunc(compress(&no_of_files, ";"));

    %let selec_line=;
    %let line_middle =;
    %let out_tabl_nam = %lowcase(&ligID)_close_contact;

    %do i = 1 %to %eval(&no_of_files);        
           
        %let file_ID = %scan(&dst_id_lst, &i, %str( ));
        %*let dst_name = &libry..%lowcase(&ligID)_&file_ID._&sufx;
        %let dst_name = %lowcase(&ligID)_&file_ID._ccname;
        %put &dst_name;
        %*put &i;
        %if &i = 1 %then
            %do;
                %let line1_1 = %str(proc sql;) create table &out_tabl_nam. as%str( ); 
                %let line1_2 = %str(select * from &libry..&dst_name. outer union corr)%str( );
            %end; 
               
       %else %if  &i > 1 and &i < %eval(&no_of_files) %then
            %do; 
                %*put &i;
                %let line_middle = &line_middle.%str(select * from &libry..&dst_name. outer union corr)%str( );
            %end; 
        %else %if &i = %eval(&no_of_files) %then
            %do;  
                %let line_last = &selec_line.%str(select * from &libry..&dst_name.)%str( )%str(; quit;);
            %end;
            
    %end; 
    %let selec_line = &line1_1.&line1_2.&line_middle.&line_last;
    &selec_line;
    %let selec_line =;

%mend close_contact_vertic_concat;
/*===============================================================================*/

/*===============================================================================*/
%macro create_datasets2();
    
    %let dst_id1 = p1 p2 p3 p4;
    %let dst_id2 = p1 p2 p3 p4 p5;
    
    %let liglst =;

    %let lengt = 10;
    
    %do k = 1 %to &lengt;
        %let lignm = %scan(&lignms., &k, %str( ));
        %let dst_id = &dst_id1;

        %if &lignm. in (DCT DCV EDF MEP) %then
            %let dst_id = &dst_id2;        
            %*put &dst_id;

        %close_contact_vertic_concat(&lignm, &dst_id, &lignm.lib);  
        %make_columns(&lignm);

        %residue_and_op_dst3(&lignm, &k);


        proc delete data = close_contact_all_temp;
        run;

        proc sort data = close_contact_all_op;
            by Res_ID OP_nm;
        run;

    %end;
%mend create_datasets2; 


/*===============================================================================*/

%macro residue_and_op_dst3(LigID, j);

    %if &j = 1 %then 
        %do;

            proc sql;
                create table close_contact_all_op as
                select Prt_res, count(Prt_res) as Res_count,
                count(Prt_res) as Normalize 
                from &LigID._close_contact2x            
                group by Prt_res;
                                 
            quit;

            proc stdize data=close_contact_all_op 
                out=close_contact_all_op method=RANGE;
                var Normalize;
            run;

            proc sql;
                alter table close_contact_all_op 
                add op_nm char (5), Res_ID num(5);

                update close_contact_all_op
                set op_nm =  "&LigID",
                    Normalize = Normalize * 1000,
                    Res_ID = input(scan(Prt_res,1,'/'), 5.);
            quit;

        %end;

   %else
        %do;  

            proc sql;
                create table close_contact_&j as
                select Prt_res, count(Prt_res) as Res_count,
                count(Prt_res) as Normalize 
                from &LigID._close_contact2x            
                group by Prt_res;
                                 
            quit;

            proc stdize data = close_contact_&j 
                out = close_contact_&j method=RANGE;
                var Normalize;
            run;

            proc sql;
                alter table close_contact_&j 
                add op_nm char (5);

                update close_contact_&j
                set op_nm =  "&LigID",
                    Normalize = Normalize * 1000;

                create table close_contact_all_temp as 
                select * 
                from close_contact_all_op
                outer union corr
                select * 
                from close_contact_&j;
            quit;


            data close_contact_all_op;
                set close_contact_all_temp;
                length Res_ID 5.;
                Res_ID = scan(Prt_res,1,'/');
            run;

            proc delete data = close_contact_&j;
            run;
        %end;

    proc sort data  =  close_contact_all_op;
        by Res_ID op_nm;
    run;

%mend residue_and_op_dst3;
/*===============================================================================*/




/*===============================================================================*/
%macro make_columns(ligID);
     /*&libry..*/
    data &ligID._close_contact2x (drop = dist_name dist);
        set &ligID._close_contact;
        where scan(dist_name,1,'_') ^= '';

        length OP_atm $3 Prt_res $7 /*Prt_atm */ opatm_prtres $15 /*opatm_prtatm $15 */ Res_ID 4.;
        *Prt_atm = scan(dist_name,2,'_');
        OP_atm = scan(scan(dist_name,1,'_'),5,'/');
        *Prt_atm = scan(scan(dist_name,2,'_'),5,'/');
        Res_ID = scan(scan(dist_name,2,'_'),4,'/');
        Prt_res = scan(scan(dist_name,2,'_'),4,'/')|| '/' || scan(scan(dist_name,2,'_'),3,'/');
        *Prt_atm = right(Prt_res)||'/'||scan(scan(dist_name,2,'_'),5,'/');
        opatm_prtres = right(OP_atm)||'-'||Prt_res;
        *opatm_prtatm = right(OP_atm)||'-'||Prt_atm;        
    run;

    proc sort data = &LigID._close_contact2x;
        by Res_ID;
    run;


%mend make_columns;

/*===============================================================================*/
/* Define some ODS template*/

/*
proc template;
    source styles.mystyle;
   source styles.htmlblue;
   source styles.statistical ;
   source styles.default;
run;
*/

%create_datasets2();

proc template ;
    define style styles.mystyle; parent = styles.HTMLBlue;
        style GraphFonts /
            'GraphValueFont' = ("<MTserif>, Times New Roman", 8pt)
            'GraphLabelFont' = ("<MTserif>, Times New Roman", 10pt, bold)
            'GraphTitleFont' = ("<MTserif>, Times New Roman", 12pt, bold)
            'GraphDataFont' = ("<MTserif>, Times New Roman", 8pt)
            'NodeLabelFont' = ("<MTserif>, Times New Roman", 10pt);
         /*style GraphDataDefault / linethickness = 2px;*/
         /*style GraphAxisLines / linethickness = 2px;*/
         /*style GraphWalls / lineThickness = 2px FrameBorder = on;*/
         /*style graphdata1 / linestyle = 1 ContrastColor = white;*/
         /*style graphdata2 / linestyle = 1 ContrastColor = black;*/
    end; 
run;
ods graphics on / reset = all;

ods pdf style= mystyle file = "&path./ALL_close_contacts_residue_1.pdf" dpi=300;
ods html style = mystyle path = "&path" gpath = "&path" dpi=300/*(url="png/")*/ 
        file = "All_close_contacts_Residues_1.htm";

ods graphics on /border=off width= 9 in height= 5 in
        ANTIALIASMAX=4862900 maxlegendarea=30 imagename = "ALL_close_contacts_residues_1" imagefmt= png;

/*===============================================================================*/

proc template;
  define statgraph gradientplot;
  dynamic _X _Y _Z _T _F _W _H;
  *mvar LEGENDTITLE "optional title for legend";
    begingraph /designwidth= _W designheight= _H; 
      *entrytitle _T; 
      layout overlay / 
        XAXISOPTS = (label="AChE residue" griddisplay=on 
                    discreteopts =(tickvaluefitpolicy=ROTATE tickvaluerotation=vertical) 
                    TICKVALUEATTRS=(size=8)) 
        YAXISOPTS = (label="OP ID" griddisplay=on);       
        heatmapparm x=_X y=_Y  colorresponse = _Z  / 
        colormodel=(lightcyan VLIGB  BIGB  VIGB STB  VIB) name = "heatmap";
        *colormodel=(whitesmoke lightcyan lightblue lightskyblue blue) name = "heatmap";
        continuouslegend "heatmap" / title="Rescaled Contact Frequency" 
            orient = horizontal location= outside halign=center 
            valign=bottom valuecounthint=10;
        
      endlayout;
      *entryfootnote _F;    
    endgraph;
  end;
run;
/*===============================================================================*/

proc sgrender data=close_contact_all_op template=gradientplot;

    dynamic _X='prt_res' _Y='OP_nm' _Z='Normalize' T='' _W= '9 in' _H = '5 in' 
        _F="";
run;




proc template;
   delete Styles.mystyle;
run;


proc sgrender data=close_contact_all_op template=gradientplot;

    dynamic _X='prt_res' _Y='OP_nm' _Z='Res_count' T='' _W= '10.5 in' _H = '5 in' 
        _F="";
run;



/*

proc sql;
    create table contact_catog as
    select Prt_res,
            (select op_nm 
            from Close_contact_all_op
            where Normalize < 200
            ) as low_contact
    from Close_contact_all_op
    group by Prt_res, low_contact;
quit;





proc format;
    value catgry 0.0 -<199.99999 = 'Low'
                  200 -< 499.99999 = 'Med_low'
                  500 -< 749.99999 = 'Med_strong'
                  750 -< 1000.999 = 'Strong';
run;

proc print data = Close_contact_all_op;
    format Normalize catgry.;
    var Prt_res Normalize op_nm;
run; 

proc sql;
    select Prt_res, Normalize format catgry., op_nm
    from Close_contact_all_op;
quit;

*/
data mltplt.contact_catogerized;
    length prt_res $10. low $50. low_med $50. med_stro $50. strong $50.;
    
    do until (last.Prt_res);
       set Close_contact_all_op;
       by Prt_res notsorted;
       if Normalize < 249.99999 then low = catx(',',low,op_nm);
       if Normalize > 250 and  Normalize < 499.99999 then low_med = catx(',',low_med,op_nm);
       if Normalize > 500 and  Normalize < 749.99999 then med_stro = catx(',',med_stro,op_nm);
       if Normalize > 750 and  Normalize < 1000.99999 then strong = catx(',',strong,op_nm);
    end;
    drop Normalize Res_count op_nm Res_ID;
run;


data mltplt.contact_catogezd_no_low;
    set mltplt.contact_catogerized;
    drop low;
    if low_med ='' and med_stro ='' and strong ='' then
        delete;
run;


proc export data=mltplt.contact_catogerized
   outfile="&path/contact_catogerized.csv"
   dbms=dlm;
   delimiter=',';
   *dbms = csv;
run;

proc export data=mltplt.contact_catogezd_no_low
   outfile="&path/contact_catogerized_no_low.csv"
   dbms=dlm;
   delimiter=',';
run;




ods html close;
ODS PDF CLOSE;
