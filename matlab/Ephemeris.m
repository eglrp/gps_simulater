function orbit = Ephemeris(sign_set)

orbit=zeros(length(sign_set.PRNmat),20);
eph=rinexe(sign_set.ephfile);
for i=1:1:length(sign_set.PRNmat)
    orbit(i,:)=eph(sign_set.PRNmat(i),2:size(eph,2));
    if sum(orbit(i,:))==0
        fprintf('error:\n�ļ��в������� %d �����ǵĹ������\n',sign_set.PRNmat(i));
    end
end