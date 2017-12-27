function eph = rinexe(ephemerisfile)
%RINEXE Reads a RINEX Navigation Message file and
%	     reformats the data into a matrix with 21
%	     rows and a column for each satellite.
%	     The matrix is stored in outputfile
%     读取32颗星在time所指定时刻的星历
%     time是一个1x6的数组，数组元素依次代表年、月、日、时、分、秒
%Typical call: rinexe('pta.96n','pta.nav')
% Revision by Wu Ling at 2017-9-30

% Units are either seconds, meters, or radians
%fide = fopen('G:\matlab\GNSS_INS组合相对导航\RINEX\brdc1980.17n');
%fide = fopen('G:\matlab\GNSS_INS组合相对导航\RINEX\brdc0060.10n');
fide = fopen(ephemerisfile);
head_lines = 0;
while 1  % We skip header
   head_lines = head_lines+1;
   line = fgetl(fide);
   answer = findstr(line,'END OF HEADER');
   if ~isempty(answer), break;	end;
end;
head_lines;
noeph = -1;
while 1
   noeph = noeph+1;
   line = fgetl(fide);
   if line == -1, break;  end
end;
noeph = noeph/8;                     % 计算一共有多少组数据
frewind(fide);
for i = 1:head_lines, line = fgetl(fide); end;

% Set aside memory for the input
svprn	 = zeros(1,32);
weekno	 = zeros(1,32);
tgd	 = zeros(1,32);
iodc	 = zeros(1,32);
toe	 = zeros(1,32);
af2	 = zeros(1,32);
af1	 = zeros(1,32);
af0	 = zeros(1,32);
aode	 = zeros(1,32);
deltan	 = zeros(1,32);
M0	 = zeros(1,32);
ecc	 = zeros(1,32);
roota	 = zeros(1,32);
toe	 = zeros(1,32);
cic	 = zeros(1,32);
crc	 = zeros(1,32);
cis	 = zeros(1,32);
crs	 = zeros(1,32);
cuc	 = zeros(1,32);
cus	 = zeros(1,32);
Omega0	 = zeros(1,32);
omega	 = zeros(1,32);
i0	 = zeros(1,32);
Omegadot = zeros(1,32);
idot	 = zeros(1,32);
svaccur = zeros(1,32);
svhealth = zeros(1,32);
L2codes = zeros(1,32);
L2flag = zeros(1,32);
tom = zeros(1,32);
IODE =zeros(1,32);
toc   =zeros(1,32);

for i = 1:noeph
   line = fgetl(fide);	  %
   temp_svprn=str2num(line(1:2));
   year = str2num(line(3:6));
   month = str2num(line(7:9));
   day = str2num(line(10:12));
   hour = str2num(line(13:15));
   minute = str2num(line(16:18));
   second = str2num(line(19:22));
   temp_time=[year month day hour minute second];     %读取指定时间对应的所有卫星星历
   if i==1
      svprn(i) = temp_svprn;
      time=temp_time;
      toc(svprn(i))=temp_time(6);
   elseif  time==temp_time
      svprn(i) = temp_svprn;
      toc(svprn(i))=temp_time(6);  
   else 
       break;
   end
   af0(svprn(i)) = str2num(line(23:41));
   af1(svprn(i)) = str2num(line(42:60));
   af2(svprn(i)) = str2num(line(61:79));
   line = fgetl(fide);	  %
   IODE(svprn(i)) = str2num(line(4:22));
   crs(svprn(i)) = str2num(line(23:41));
   deltan(svprn(i)) = str2num(line(42:60));
   M0(svprn(i)) = str2num(line(61:79));
   line = fgetl(fide);	  %
   cuc(svprn(i)) = str2num(line(4:22));
   ecc(svprn(i)) = str2num(line(23:41));
   cus(svprn(i)) = str2num(line(42:60));
   roota(svprn(i)) = str2num(line(61:79));
   line=fgetl(fide);
   toe(svprn(i)) = str2num(line(4:22));
   cic(svprn(i)) = str2num(line(23:41));
   Omega0(svprn(i)) = str2num(line(42:60));
   cis(svprn(i)) = str2num(line(61:79));
   line = fgetl(fide);	    %
   i0(svprn(i)) =  str2num(line(4:22));
   crc(svprn(i)) = str2num(line(23:41));
   omega(svprn(i)) = str2num(line(42:60));
   Omegadot(svprn(i)) = str2num(line(61:79));
   line = fgetl(fide);	    %
   idot(svprn(i)) = str2num(line(4:22));
   L2codes(svprn(i)) = str2num(line(23:41));
   weekno(svprn(i)) = str2num(line(42:60));
   L2flag(svprn(i)) = str2num(line(61:79));
   line = fgetl(fide);	    %
   svaccur(svprn(i)) = str2num(line(4:22));
   svhealth(svprn(i)) = str2num(line(23:41));
   tgd(svprn(i)) = str2num(line(42:60));
   iodc(svprn(i)) = str2num(line(61:79));
   line = fgetl(fide);	    %
   tom(svprn(i)) = str2num(line(4:22));
   i;
   end
status = fclose(fide);

%  Description of variable eph.
eph =  zeros(32,21);
eph(:,1) = linspace(1,32,32);    % 32颗卫星编号
eph(:,2)  = crs;
eph(:,3)  = deltan;
eph(:,4)  = M0;
eph(:,5)  = cuc;
eph(:,6)  = ecc;
eph(:,7)  = cus;
eph(:,8)  = roota;
eph(:,9)  = toe;
eph(:,10)  = cic;
eph(:,11) = Omega0;
eph(:,12) = cis;
eph(:,13) = i0;
eph(:,14) = crc;
eph(:,15) = omega;
eph(:,16) = Omegadot;
eph(:,17) = idot;
eph(:,18) = af0;
eph(:,19) = af1;
eph(:,20) = af2;
eph(:,21) = toc;
%fidu = fopen(outputfile,'w');
%count = fwrite(fidu,[eph],'double');
fclose all;
%%%%%%%%% end rinexe.m %%%%%%%%%
