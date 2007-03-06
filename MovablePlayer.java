/**
 * A <code>MovablePlayer</code> is an object which is used to play the game,
 * it is possible to read and write data from this player
 *
 * @author <a href="mailto:m.w.manders@student.tue.nl">Maarten Manders</a>
 */
public class MovablePlayer extends ReadablePlayer
{
    /**
     * Creates a new <code>MovablePlayer</code> object, and
     * assigns a given id to it.
     *
     * @param id The id to add to this <code>MovablePlayer</code>.
     * @param nickname The nickname of this player.
     */
    public MovablePlayer(int id, String nickname)
    {
        super(id, nickname);
    }

    /**
     * Sets the horizontal position of this player.
     *
     * @param x_pos The new horizontal position
     */
    public void setXPos(int x_pos)
    {
        super.setXPos(x_pos);
    }

    /**
     * Sets the vertical position of this player.
     *
     * @param y_pos The new vertical position
     */
    public void setYPos(int y_pos)
    {
        super.setYPos(y_pos);
    }
}