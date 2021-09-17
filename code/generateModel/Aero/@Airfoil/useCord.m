function AF = useCord(foilcords)
    
    AF = Airfoil;
    AF.UpperX = foilcords(end:-1:1,1);
    AF.UpperY = foilcords(end:-1:1,2);
    AF.LowerX = foilcords(:,3);
    AF.LowerY = foilcords(:,4);
    AF.Name   = 'Airfoil';
    
end